---
name: hit-triage
description: >
  This skill should be used when the user asks to "triage these hits", "cluster these compounds",
  "how many series are in this hit list", "find the common core", "pick compounds to follow up",
  "which analogues should I order", "group by scaffold", "find neighbours of this hit", or runs
  /hit-triage. Covers clustering into chemotypes, scaffold analysis, MCS/series cores, similarity
  searching and representative selection.
metadata:
  version: "0.1.0"
---

# Triaging a hit list

A screen returns compounds; a chemist needs *series*. The questions that matter are how many distinct
chemotypes are really there, which ones have support from more than one member, and which compounds
to order next. Singletons and clusters are judged completely differently.

## Procedure

1. **Clean and deduplicate first.** `standardize_library`. Hit lists from screening decks routinely
   contain salt forms and repeats of the same compound; counting those as separate hits inflates
   every number that follows.
2. **Cluster into chemotypes.** `cluster_library`, default cutoff 0.35 distance (0.65 Tanimoto on
   ECFP4). Report cluster count, size distribution and singleton count. Then state the thing that
   matters: **a cluster of several related actives is evidence of SAR; a singleton is a compound that
   has not reproduced.** Singletons are where frequent hitters and assay artefacts concentrate.
3. **Check the scaffold view too.** `scaffold_analysis` groups by Bemis-Murcko scaffold, which is a
   stricter notion than fingerprint similarity. Run it with `generic=false` first, then
   `generic=true` to merge bioisosteric variants of one core. Where the two views disagree, say so -
   it usually means a series is held together by a side chain rather than a core.
4. **Find the core of each real series.** `series_core` on the members of a cluster returns the MCS.
   Use `threshold=0.8` when one member is structurally off. Then depict the series with
   `align_to` set to that core, so the R-groups line up and the SAR is readable at a glance. This is
   the single most useful artefact to hand a chemist.
   MCS over a whole heterogeneous deck is meaningless - it returns a trivial fragment. Only run it
   within a series.
5. **Expand the interesting hits.** `similarity_search` against the full deck or a vendor file pulls
   the near neighbours that were screened but did not hit. Neighbours that were tested and were
   *inactive* are as informative as the hit: a cluster with a flat SAR across close analogues is
   suspicious, while a sharp, sensible SAR is the best evidence a series is real.
6. **Choose what to follow up.** `cluster_library` with `pick_per_cluster` gives a diverse
   representative set (MaxMin within each cluster). Rank clusters for the user by: number of members,
   SAR consistency, property profile of the best member, absence of alerts, and synthetic
   tractability. Say which criterion drove each ranking.
7. **Screen the triaged set for properties and alerts.** Hand off to `/screen-library` and
   `/alerts-review` rather than repeating that work here.

## Reporting

Deliver: cluster count and sizes, the series worth pursuing with their cores drawn and aligned, the
singletons listed separately with the reason they are riskier, a recommended follow-up set with the
basis for each pick, and the files written.

Be explicit about what clustering does and does not mean. Tanimoto similarity on a fingerprint is a
structural bookkeeping measure, not a biological one: two compounds at 0.9 similarity can differ
tenfold in potency (an activity cliff), and two at 0.3 can hit the same site. The cutoff is a
convention, so when a decision turns on it, try a second cutoff and report whether the picture holds.

Never rank hits by potency alone, and never call a cluster "the best series" on structure alone -
that judgment needs assay quality, selectivity, orthogonal confirmation and the chemist's knowledge
of the target. Lay out the evidence and say what you would want to see next.
