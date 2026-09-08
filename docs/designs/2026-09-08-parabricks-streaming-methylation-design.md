# Low-latency methylation annotation in Parabricks `fq2bam_meth`

**Status:** approved design; implementation not started
**Decision date:** 2026-09-08

## Decision

Add an opt-in `fq2bam_meth` mode, proposed as `--emit-methylation-xm`, that
adds an `XM:Z` tag to **every mapped SAM/BAM alignment record**. The tag is
computed as soon as that record has its final CIGAR, aligned coordinate, and
bisulfite-conversion orientation, before BAM serialization and before the
sort, duplicate-marking, and optional BQSR stages.

This is an annotation stage, not a site-level methylation aggregator. It must
be CIGAR-aware and run on the GPU whenever the alignment record is
GPU-resident. The tag remains attached to its record as subsequent stages
reorder records or change flags/qualities.

## Why this boundary

`fq2bam_meth` is documented as a GPU implementation of bwa-meth using BWA-MEM
as its alignment backend, followed by coordinate sorting, optional duplicate
marking, and optional BQSR. Its `--align-only` mode emits BAM after BWA-MEM,
before sorting and duplicate marking. Upstream BWA generates a CIGAR only once
alignment endpoints are known.

Therefore:

* A seed or provisional chain is too early: local extension, split alignment,
  and primary/secondary selection can change the coordinate or CIGAR.
* A final alignment record is sufficient: it contains the query sequence,
  final coordinate, CIGAR, flags, mapping quality, and conversion orientation.
* Sorting and duplicate marking are too late for the desired latency and do
  not alter the aligned bases or CIGAR. BQSR alters base qualities, but not
  the sequence, CIGAR, or reference context.

The main uncertainty is implementation placement, not the contract: public
documentation does not reveal whether Parabricks materializes final CIGARs in
a device-resident record, a host record, or in the serializer. The Parabricks
source owners must locate that materialization boundary.

## Scope

### In scope

* `XM:Z` annotation for all mapped alignment records: primary, secondary,
  supplementary, duplicate-marked, and low-MapQ records.
* CpG, CHG, CHH, unknown-context, and non-informative per-query-base symbols.
* A GPU annotation stage between finalized alignment and the output/sort queue.
* Exact parity tests against a CPU reference annotator on the final record.
* Preservation of `XM` through Parabricks sorting, duplicate marking, and
  optional BQSR.

### Out of scope

* GPU site-level aggregation, methylation percentages, or deduplication of
  methylation calls. Those need downstream policy and are not intrinsic to an
  individual alignment.
* Replacing MethylDackel's extract/mbias workflows.
* Revising an `XM` tag after BQSR.
* New long-read behavior. If the selected Parabricks configuration falls back
  to CPU for a read, it must use the same annotator semantics or fail the
  opt-in mode; it must never silently omit `XM` for a mapped record.

## Record and tag contract

### Command-line interface

Proposed opt-in flag:

```text
--emit-methylation-xm
```

It requires no additional reference argument. `fq2bam_meth` already receives
the original reference through `--ref` and locates the converted
`.bwameth.c2t` reference for alignment. The annotator must inspect the
**original, unconverted** reference sequence.

For a mapped record, append exactly one `XM:Z:<calls>` tag. `<calls>` has one
ASCII character per stored `SEQ` base (the query length after hard clipping).
The stage must not overwrite a pre-existing `XM` tag; that condition is an
internal pipeline error because this FASTQ-to-BAM workflow owns the output
record.

The final record must also carry, or make available to the annotator, the
bwa-meth conversion orientation. Emit/retain `YD:Z:f` for the C-to-T/original
top orientation and `YD:Z:r` for the G-to-A/original-bottom orientation. This
is the convention already recognized by MethylDackel.

Unmapped records have no reference coordinate and receive no `XM` tag.

### Symbol semantics

| Symbol | Meaning |
| --- | --- |
| `Z` / `z` | methylated / unmethylated CpG observation |
| `X` / `x` | methylated / unmethylated CHG observation |
| `H` / `h` | methylated / unmethylated CHH observation |
| `U` / `u` | methylated / unmethylated cytosine whose context cannot be determined |
| `.` | no methylation observation for this query base |

Uppercase/lowercase represent the observed methylated/unmethylated base in
the conversion orientation: C/T for `YD:f`, G/A for `YD:r`. A mismatch,
ambiguous query base, insertion, soft-clipped base, or other non-callable
query position is `.`. At a contig boundary, or where the needed neighbouring
reference bases are ambiguous, use `U`/`u` when the observed base supports a
cytosine call but the CpG/CHG/CHH class cannot be determined.

`XM` is intentionally **not quality-filtered**. Every mapped record receives
the full reference-context annotation. Consumers apply base-quality, MapQ,
duplicate, primary/secondary, and supplementary-alignment policy using the
record's final `QUAL` and flags. This keeps `XM` stable if optional BQSR later
changes qualities and avoids dropping information at annotation time.

## Annotation algorithm

For each finalized mapped alignment record:

1. Validate that its CIGAR consumes the stored query length and its reference
   span lies on the selected original-reference contig. Validate conversion
   orientation.
2. Allocate a query-length call buffer and initialize it to `.`.
3. Walk the final CIGAR while maintaining query and reference offsets.
   `M`, `=`, and `X` advance both and map a query base to a reference base.
   `I` and `S` advance only the query and remain `.`. `D` and `N` advance only
   the reference. `H` and `P` advance neither.
4. At every query-to-reference position, inspect the original-reference base
   and the required flanking bases. Classify CpG, CHG, CHH, or unknown context,
   then emit the symbol warranted by the oriented observed base.
5. Attach `XM:Z` and `YD:Z` to the same finalized record and submit the record
   to the existing output/sort queue.

The algorithm must use 64-bit-safe reference coordinates internally and must
not reconstruct a reference-aligned read sequence. That avoids the insertion,
deletion, soft-clip, and hard-clip errors fixed in MethylDackel's current
per-read mapping logic.

## Pipeline placement and data flow

```text
FASTQ + original reference + converted bwa-meth index
    -> candidate finding / chaining / local extension
    -> final alignment record (CIGAR, position, flags, MapQ, YD)
    -> GPU XM annotator
    -> BAM record encoder / existing asynchronous sort queue
    -> optional duplicate marking -> optional BQSR -> BAM or CRAM
```

The annotator must consume the existing internal final-alignment batch rather
than parse serialized BAM. It should reuse the original-reference
representation already available to the alignment pipeline, or add a read-only
original-reference accessor/cache. Its output allocation budget is one byte
per stored query base plus SAM auxiliary-tag overhead; no global
coordinate-order state is required.

The exact source-level hook is to be identified by the Parabricks owners:

1. Find the code that turns a selected extended alignment into final CIGAR,
   coordinates, flags, MapQ, and bwa-meth orientation.
2. Confirm whether that batch is device-resident. Add the annotation kernel
   there; otherwise add a vectorized CPU implementation at that boundary.
3. Ensure the serializer accepts the appended variable-length `XM` and `YD`
   tags without changing ordering or dropping tags in later stages.

## Filtering and BQSR policy

Annotation has no MapQ, base-quality, duplicate, or primary-only filter.
Those are downstream interpretation policies:

* `XM` is present for every mapped alignment, including secondary and
  supplementary records; a multi-mapping read may therefore have one distinct
  tag on each alignment record.
* Duplicate marking may set a flag after annotation. Aggregators exclude or
  include that record according to their normal duplicate policy.
* BQSR may change `QUAL`. Consumers that filter methylation observations use
  the final `QUAL` with the query index of the corresponding `XM` character.
* The default consumer policy should remain explicit: for example, primary,
  non-duplicate alignments with MapQ >= 10 and a chosen base-quality cutoff.
  The tag itself makes no such decision.

This differs deliberately from MethylDackel `perRead`, whose current emitted
string applies its requested base-quality cutoff during construction. The
Parabricks tag is a stable per-base annotation; MethylDackel remains the CPU
reference for context/orientation rules and may later gain an optional
tag-consuming mode.

## Failure handling and observability

Fail the opt-in invocation rather than write a misleading tag if a mapped
record has a malformed CIGAR, an unavailable/mismatched original-reference
contig, an out-of-range alignment span, missing conversion orientation, or an
auxiliary-tag encoding failure. Contig-edge and ambiguous-neighbour contexts
are valid biological/reference cases and produce `U`/`u`, not errors.

Expose counters in the existing tool metrics/logging surface:

* mapped records annotated;
* unmapped records passed through;
* query bases marked CpG, CHG, CHH, unknown, and `.`;
* CPU-fallback records annotated;
* failures by validation class.

At startup, validate original-reference contig names and lengths against the
alignment reference dictionary. Where reference MD5s are available, validate
them as well. This prevents context calls against the converted or wrong
reference.

## Validation and acceptance criteria

### Correctness

Build a small, independently implemented CPU oracle that takes a final SAM
record, its original reference, and its `YD` value and emits `XM`. Compare the
GPU result byte-for-byte per record. The oracle must not use output BAM order.

The conformance corpus must cover:

* `YD:f` and `YD:r` and all `Z/z`, `X/x`, `H/h`, `U/u`, and `.` symbols;
* insertions, deletions, skipped regions, soft clips, hard clips, padding, and
  contig boundaries;
* primary, secondary, supplementary, duplicate-marked, low-MapQ, and unmapped
  records;
* a representative random sample of scrubbed, actual read sequences in
  addition to synthetic cases;
* BQSR enabled and disabled (identical `XM` for the same pre-BQSR record);
* duplicate marking enabled and disabled (identical `XM` for matching records);
* the CPU fallback path, if it is allowed by the selected read-length mode.

The existing self-contained MethylDackel per-read fixtures are suitable seeds
for the synthetic CIGAR/context cases. The Parabricks conformance test must
use only synthetic or scrubbed records and references.

### Performance

Benchmark the tagged path against the same `fq2bam_meth` command without the
flag on representative WGBS data and GPU configurations. The initial release
gate is no more than 5% end-to-end throughput regression and no unbounded
device-memory growth. Report annotation-stage time, bytes of tag output, and
host/device transfer changes separately. If this gate is missed, retain the
same contract but optimize the kernel/data layout; do not move annotation to a
post-BAM MethylDackel pass as a silent substitute.

## Rollout

1. **Source spike:** identify final-CIGAR materialization, reference access,
   conversion-orientation representation, tag serializer support, and the
   CPU-fallback path. Confirm the proposed boundary with a one-read trace.
2. **Feature implementation:** add the opt-in flag, record annotator, output
   tag encoder, counters, and CPU oracle/conformance tests.
3. **Parity and performance gate:** run the corpus and benchmark matrix; fix
   discrepancies before exposing the flag beyond preview use.
4. **Preview release:** document tag semantics and downstream filtering policy.
   Preserve existing output exactly when the flag is absent.
5. **Consumer follow-up (optional):** add a MethylDackel mode that reads `XM`
   while applying user-selected quality/flag filters, avoiding recomputation.

## Alternatives rejected

* **Annotate from seeds or provisional chains:** rejected because no stable
  CIGAR/coordinate or primary-alignment decision exists.
* **Run MethylDackel only after sorted BAM output:** correct but loses pipeline
  overlap, requires serialized BAM and often an index, and adds storage/I/O.
* **Annotate only primary, high-quality reads:** rejected because it destroys
  data and bakes downstream filtering policy into a lossless record tag.
* **Apply a base-quality cutoff while writing `XM`:** rejected because optional
  BQSR can change that cutoff decision. Store the full annotation and filter
  at consumption.

## References

* NVIDIA, [fq2bam_meth documentation](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam_meth): GPU bwa-meth-compatible BWA-MEM backend, sorting, optional duplicate marking/BQSR, and `--align-only` behavior.
* NVIDIA, [fq2bam_meth release note](https://docs.nvidia.com/clara/parabricks/4.5.0/whatsnew/newtools_v4.3.0-1.html): compatible pre/post-processing around BWA-MEM and shared accelerated alignment code.
* Heng Li, [BWA source](https://github.com/lh3/bwa/blob/master/bwa.c): CIGAR generation after alignment endpoints are known.
