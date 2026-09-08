# Upstreaming decomposition

How the fork's work can be offered to `jltsiren/gcsa2` and `vgteam/vg` as
separately reviewable, separately acceptable pieces, and what stands in the way
of that today.

This is a companion to [EXTERNAL_MEMORY_CONSTRUCTION.md](EXTERNAL_MEMORY_CONSTRUCTION.md),
which is the status of record for what the external route does. This document
is about packaging, not behavior.

## Why the question is not trivial

Every optimization here was gated on byte-identical output. That makes each one
individually defensible, and it invites the assumption that they are therefore
individually *takeable*. They are not, and for two different reasons in the two
repositories.

In `vg`, the changes are separable but were never separated: five unrelated
performance fixes, a GCSA2 integration, and an mpmap research feature all sit on
one branch, in one line of history. Nothing is entangled; nothing is packaged.

In `gcsa2`, the situation is the reverse. The work is packaged as a coherent
feature but is genuinely inseparable: 74 of 86 commits introduce new files that
later commits build on, and the three shared primitives that the legacy route
also uses were modified inside the very first commit.

The practical consequence is that "pick and choose" means something different in
each repository, and a single merge strategy cannot serve both.

## Measured divergence surface

Measured against `e90ba1a` (upstream gcsa2 v1.4.0) and `7272a05e7` (the vg merge
base). Line counts are `git diff --numstat`; classification of new-versus-
modified is by `git cat-file -e` against the base tree.

### gcsa2

86 commits, no merge commits, +26,999 / −311 across 57 files. 37 of those files
are new; 20 existed upstream.

The deletion count is the number that matters. 311 deleted lines against 27,000
added means the fork is overwhelmingly additive, and the deletions concentrate
almost entirely in three files:

| File | Added | Deleted | Share of all deletions |
|---|---|---|---|
| `src/path_graph.cpp` | 3,151 | 206 | 66% |
| `src/gcsa.cpp` | 1,278 | 41 | 13% |
| `include/gcsa/internal.h` | 149 | 39 | 13% |
| everything else (17 files) | 2,722 | 25 | 8% |

The dispatch seam is narrow. `ConstructionParameters::externalMemory()` is
consulted at 17 sites in the whole library — 8 in `gcsa.cpp`, 4 in
`build_gcsa.cpp`, and one each in `lcp.cpp`, `files.cpp`, and `algorithms.cpp`,
plus 2 in `support.h` itself. Header coupling is narrower still: exactly two
upstream headers gained an include of a fork header (`path_graph.h` pulls in
`workspace.h`; `algorithms.h` pulls in `external_sort.h`).

### vg

96 commits, +8,419 / −201 across 56 files, excluding submodules. Three
unrelated bodies of work share the branch:

- **GCSA2 integration** (~7 commits): `index_registry.cpp`, `index_main.cpp`,
  `autoindex_main.cpp`, `gcsa_workspace.cpp`, `gcsa_helper.cpp`,
  `gcsa_worker_main.cpp`. This is feature plumbing and cannot precede the gcsa2
  work.
- **mpmap splice research** (~16 commits on `multipath_mapper.cpp` alone):
  `--trace-splice-search`, MMP seeding, `--splice-denovo`, junction
  canonicalization. This is not output-neutral and is not performance work. It
  does not belong in any conversation about these optimizations.
- **Five standalone performance fixes**, each confined to one subsystem in one
  commit.

## Modularity

The intrusion into upstream code falls into three tiers, and a maintainer's risk
is entirely determined by which tier a change sits in.

### Tier 1 — dispatch seam (low risk)

New route branched around untouched upstream code. This describes the 17
`externalMemory()` sites, all 37 new gcsa2 files, and the `path_graph_external.h`
interface, which exposes the entire external path-graph subsystem through 205
lines and 14 free functions. `src/external_join.cpp` is 5,109 lines but
holds 19 types inside an anonymous namespace; its public surface is small and its
boundary is real. A maintainer can review this tier by reading the seam, not the
implementation.

### Tier 2 — additive change to shared structures (moderate risk)

`GCSA::GCSA(InputGraph&, ...)` grew from 294 lines upstream to 555, with 8
external-route branches interleaved through the body. Nothing upstream was
removed; the legacy path still reads top to bottom. But the function is now
twice as long as it was and no longer has a single subject, and every future
upstream change to construction has to be made around eight conditionals. The
fork already began undoing this — `ce56be9` lifted external reporting out into
`reportPruneMergeStats`, `reportJoinStats`, `reportFinalMergeStats`, and
`reportFinalEventStats` — but the branching itself remains inline.

`ConstructionParameters` gained 18 setter/getter pairs and `build_gcsa` gained
about 20 CLI flags, all flat on a class upstream keeps deliberately small.

### Tier 3 — in-place rewrite of code the legacy route runs (the real risk)

This is the tier that undermines the "legacy route is preserved byte-for-byte"
claim, which is true of the *output* but not of the *code path*.

Three shared structures were rewritten rather than branched around:

- **`SameFromSet`** (`path_graph.cpp`) was replaced entirely. Upstream's
  `std::vector<node_type> nodes, buffer` members and its `fromNodes` /
  `operator()` / `select` implementation are all gone.
- **`PathGraphBuilder`** lost its `std::vector<WriteBuffer<PathNode>> path_files`
  and `rank_files` members in favour of a new `PathGraphOutputCache`, which caps
  concurrently open writer pairs at a default of 16. The constructor gained three
  defaulted parameters, so legacy call sites still compile unchanged — but they
  now run through bounded-fanout file handling that upstream did not have.
- **`ReadBuffer`** (`internal.h`) moved from `std::ifstream` to a raw POSIX
  descriptor, to allow `pread()` from the background reader thread and
  `posix_fadvise()` retirement of consumed ranges. Its sizing changed from
  `READ_BUFFER_SIZE = MEGABYTE` **elements** to `DEFAULT_BUFFER_BYTES = MEGABYTE`
  **bytes**.

That last one deserves stating precisely, because it is the clearest case of
"output-neutral" being read as more than it means. `PathNode` is 24 bytes
(`node_type from, to` plus `size_type fields`). Upstream reserved 1,048,576
elements, or about 25 MB per `PathNode` stream. The fork reserves 1 MB, or 43,690
elements — a 24-fold reduction in buffered elements. `MergedGraphReader::init`
defaults `buffer_bytes` to `DEFAULT_BUFFER_BYTES`, and the merged-graph scan runs
on **both** routes. The legacy in-memory route therefore has a materially
different memory and syscall profile than upstream v1.4.0, even though it emits
identical bytes.

The change is deliberate and the code says why: a `PathNode` reader was reserving
many times what a byte reader reserved, multiplied per file. The objection is not
that it is wrong. It is that a maintainer evaluating "does this touch the route I
already ship?" gets the answer "no" from the commit messages and "yes" from the
code, and the only evidence that the legacy route is unharmed is that its output
did not change — which was never the thing at risk.

Compounding this, `internal.h` was modified in `2290e53`, the *first* fork
commit, bundled with the durable workspace and external path sorting. There is no
commit that isolates the shared-primitive change. `src/path_graph.cpp` is touched
by 20 separate commits spanning the whole history.

## Simplicity

Four observations, in descending order of how much they would cost a maintainer.

**The construction constructor has no single subject.** 555 lines and 8
route conditionals. The two routes share a preamble, a doubling loop, and a
finalization, but diverge inside each. Splitting the body into three named
phases, each with a legacy and an external implementation selected once, would
let a reader follow either route without reading the other.

**Budget plumbing has sprawled into the public header.**
`path_graph_external.h` exports `pathMergeInputBudget` (two overloads),
`pathMergeInputPairs` (two overloads), `pathMergeOutputPairs`,
`pathMergeCeilingBudget`, `pathMergeInputCacheBudget`,
`mergeAdmissibleBlockSize`, `externalPathGraphSortMinimumBudget`,
`externalPathJoinMinimumBudget`, and `externalPathGraphShardPeakBytes`. Nine distinct budget entry points across eleven overloads, for one subsystem. These are derived quantities of a single
resource plan; the plan should be the exported type and the arithmetic should be
its methods.

**Eighteen flat knobs on `ConstructionParameters`.** The fork already
distinguishes semantic from operational settings (`constructionSemanticSettings`
versus `constructionOperationalSettings` in `gcsa.cpp`), which is the right
distinction and is load-bearing for workspace resume. That distinction is not
reflected in the type. One `ExternalMemoryOptions` member would add a single
field to upstream's class instead of 36 accessors, and would make the
semantic/operational split visible where it is decided rather than only where it
is consumed.

**Large translation units, but along real seams.** `external_join.cpp` at 5,109
lines, `final_events.cpp` at 2,425, `path_graph.cpp` at 4,172. `external_join.cpp`
already groups into four clusters that could each be a file: record and run I/O
(`JoinRecord`, `JoinRun`, `JoinFileWriter`, `JoinFileReader`, the
`JoinGroupSidecar` types), partition planning (`JoinPartition`, `JoinKeySample`,
`JoinRadixPack`, `JoinRadixPlan`), the sorter (`ExternalJoinSorter`,
`PathShardReader`, `JoinHeapComparator`), and worker coordination
(`ExternalJoinWorkerTask`, `ExternalJoinWorkerResult`, `ActiveJoinWorker`,
`ConcurrentPathBudgets`). This is the least urgent item: the file is large but
encapsulated, and splitting it is cosmetic next to the tier-3 issues.

## Decomposition plan

### vg: five independent pull requests, available now

Each of these is already a single atomic commit touching one subsystem. They can
be cherry-picked onto a fresh branch off `vgteam/vg` master in any order, with
one exception noted below.

| # | Commit | Files | Depends on | Notes |
|---|---|---|---|---|
| V1 | `409c30a77` bulk path deletion | `subcommand/prune_main.cpp` (+14/−4) | nothing | Calls `MutablePathHandleGraph::destroy_paths`, which already exists in upstream libhandlegraph and whose default implementation is the loop being replaced. Cannot break any graph type. |
| V2 | `0d6a46873` transcript copy elimination | `transcriptome.{cpp,hpp}` (+35/−5) | nothing | Restores implicit moves suppressed by user-declared destructors; removes two by-value range-fors over protobuf `Mapping`s. |
| V3 | `93f6f6271` streaming source-sink discovery | `source_sink_overlay.cpp` (+115/−42), unit test (+29) | nothing | Ships with its own unit test. |
| V4 | `b47de4db9` batch edge deletion | `algorithms/prune.cpp` (+41/−15), `t/38_vg_prune.t` (+13/−2), libbdsg bump | **libbdsg PR** | Needs `destroy_edges_bulk`, which is a fork addition to libbdsg (`b07563b`). A libbdsg PR must land first. |
| V5 | `bbf264574` parallel `PhaseUnfolder` | `phase_unfolder.{cpp,hpp}` (+108/−14) | nothing | **Hold.** See below. |

Offer V1 and V2 first. V1 is a three-line change against an existing upstream
API and is the strongest opening: it establishes the byte-identity discipline
with a change a reviewer can verify by inspection.

V5 should not be offered yet. Its correctness argument is that `NodeMapping::insert`
is a sequential counter, so workers unfold into local id spaces and a serial
re-base restores global ids in component order. That was validated byte-identical
at 1, 16, and 48 threads on an eight-component fixture, but has never been run at
scale — the chr21 A/B was dropped before producing a number. An eight-component
fixture is not evidence about a pangenome's component count or size distribution,
and the re-base is exactly the step whose cost grows with both. Offering it
without a scale gate asks a maintainer to accept a threading change on a fixture
that cannot exhibit the failure mode.

**Prerequisite for all five: split the branch.** `rna-copy-elimination` currently
carries GCSA2 plumbing, mpmap splice research, and these five fixes on one line
of history. No maintainer can take a piece of it. Cut five topic branches off
`vgteam/vg` master and cherry-pick; the commits are already atomic, so this is
mechanical.

### gcsa2: a stack, not a menu

The external route cannot be offered as independent optimizations, because it is
not a set of optimizations — it is one feature with a long tail of tuning. What
it *can* be is a stack of layers, each individually reviewable and each
individually revertible at an interface boundary. Proposed cut, following the
existing file structure:

- **G0 — shared primitives.** `internal.h` / `internal.cpp` (`ReadBuffer` on a
  POSIX descriptor, byte-denominated buffer sizing, `BufferWindow`),
  `resources.{h,cpp}`, `utils.{h,cpp}`, `disk_array.{h,cpp}`. **This is the only
  layer that changes code the legacy route runs**, and it must be cut as its own
  commit with its own measurement of the legacy route's memory and I/O profile —
  not just its output. Today it is buried in `2290e53`.
- **G1 — durable workspace.** `workspace.{h,cpp}`, `checkpoint.{h,cpp}`. No
  construction behavior; an artifact store with a commit/recovery protocol.
- **G2 — external sort and preprocessing.** `external_sort.*`,
  `external_preprocessing.*`, `path_sort_run.*`.
- **G3 — external join.** `external_join.cpp`, `path_graph_external.h`.
- **G4 — final events and streaming LCP.** `final_events.*`, the `lcp.cpp`
  additions, the `algorithms.cpp` external verifier.
- **G5 — framed compression.** `compressed_block.{h,cpp}`. Genuinely optional:
  the format is version-2 framed with a byte-identical version-1 logical payload,
  so a maintainer can take G0–G4 and decline G5 without touching semantics.
- **G6 — integration.** The `gcsa.cpp`, `build_gcsa.cpp`, `path_graph.cpp`,
  `files.cpp`, `support.h` changes and the CLI.
- **G7 — post-integration tuning.** The roughly 20 most recent commits
  (`90a28f2`, `ade486a`, `9c71587`, `5466d07`, `216f74b`, `2b9fc3f`, and so on).
  Each is a measured improvement to the external route and each is meaningless
  without it, but they are individually revertible and can be offered as a
  follow-up series rather than folded into G6.

G5 is the only layer a maintainer can genuinely decline while keeping the rest.
G0 is the only layer that carries risk to code they already ship. Those two facts
should lead the submission, because they are the two a reviewer most needs and
would otherwise have to derive from 27,000 lines.

Two structural changes would materially improve the odds on G6, and both are
worth doing before offering anything:

1. Collapse the 18 flat `ConstructionParameters` knobs into one
   `ExternalMemoryOptions` member, split internally along the existing
   semantic/operational boundary.
2. Split the 555-line constructor into named phases so each route is readable
   on its own.

## What this analysis does not establish

The performance figures cited here are carried forward from the session record
and were not re-measured: the 10.1x peak-memory reduction on the joint
chr20+21+22 index (64.5 GiB external against 652.3 GiB legacy, byte-identical
output), the 4.1x–27.3x range for bulk path deletion on a synthetic depth sweep,
and the 8.85% peak-RSS reduction for transcript copy elimination (147.17 → 134.15
GiB, against a same-binary noise floor of 0.027%). This document takes them as
reported.

Nothing here has been compiled or tested. The claims are structural — commit
boundaries, file classification, call-site counts, header includes, type sizes —
and every one is reproducible from the git history with the commands in this
document. No claim is made that any proposed re-cut preserves byte-identity;
re-cutting G0 in particular changes which commit contains which hunk, and the
byte-identity gates would have to be re-run against the re-cut stack.

The estimate that G5 is cleanly declinable rests on the documented invariant that
a framed join run contains the byte-identical version-1 logical header, records,
and footer. That invariant was read from the design document, not verified by
building with compression disabled.
