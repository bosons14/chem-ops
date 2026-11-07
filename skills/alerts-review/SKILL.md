---
name: alerts-review
description: >
  This skill should be used when the user asks "is this a PAINS", "check for structural alerts",
  "are these frequent hitters", "any reactive groups", "should I throw this hit out", "flag
  toxicophores", "run Brenk/NIH filters", or runs /alerts-review. Covers PAINS, Brenk, NIH, ZINC and
  ChEMBL house filters, and how to interpret an alert rather than obey it.
metadata:
  version: "0.1.0"
---

# Reviewing structural alerts

An alert says a substructure has caused trouble before. It does not say this compound is a false
positive. The literature on PAINS is routinely over-applied: the filters were derived from a
particular assay technology, and treating them as a blacklist discards real chemical matter.

Your job is to tell the chemist *which* alert fired, *why* that substructure is on the list, and what
experiment would settle it - not to recommend deletion.

## Procedure

1. **Run the catalogs.** `structural_alerts` with the compounds, `image_path` set so the flagged
   atoms are highlighted. Default catalogs are PAINS, Brenk and NIH; call `chem_tool_status` for the
   full list when the user wants vendor-style filters (ZINC, Glaxo, BMS, Dundee, SureChEMBL).
2. **Separate the three kinds of alert**, because they imply different actions:
   - **Assay interference (PAINS)** - the compound may hit in many unrelated screens. Action:
     confirm in an orthogonal assay format, check for a dose-response with a sensible Hill slope,
     test for aggregation with detergent, and look for a real SAR in the series.
   - **Reactivity (acyl halides, epoxides, Michael acceptors, aldehydes)** - may react with the
     protein or with assay components. Action: decide whether reactivity is *intended*. In a covalent
     programme an acrylamide warhead is the design, not a flaw - say so rather than flagging it as a
     problem.
   - **Metabolic or toxicophoric liability (nitro, aniline, thiophene, hydrazine, long aliphatic
     chains)** - a downstream development risk, not a reason to drop a hit at triage.
   Read `references/alert-families.md` for what each common alert family means mechanistically.
3. **Weigh the evidence before advising anything.** State the case on both sides:
   - Evidence the alert matters: the compound hits in multiple unrelated screens, flat SAR across
     close analogues, steep or shallow dose-response, activity lost with detergent or with fresh DTT,
     potency that does not track structure.
   - Evidence it does not: clean orthogonal readout, a crystal or cryo-EM structure, sensible SAR
     across a series, selectivity against related targets, reproducible potency from fresh solid.
4. **Give the counterexamples honestly.** Many marketed drugs carry alerts. Screening this plugin's
   own example file with alerts set to reject discards gefitinib and erlotinib. Aspirin trips Brenk's
   `phenol_ester`. Curcumin and quercetin trip PAINS and genuinely are promiscuous. The catalogs
   cannot tell those cases apart; the chemist and the orthogonal assay can.
5. **Recommend an experiment, not a deletion.** The useful output of this skill is the next
   experiment that would resolve the ambiguity.

## Reporting

Per compound: which alerts fired, from which catalog, on which atoms, what that substructure is
known to do, and the experiment that would settle it. Across a set: the alert histogram, since one
substructure dominating a hit list is itself a finding about the deck or the assay.

Do not tell the user to discard compounds on alerts alone. If they ask for a hard filter, apply it,
but report what it removed and name any well-known drugs or likely-genuine chemotypes in that
discard pile.

State plainly that these catalogs flag *known* problem substructures: a clean result means nothing
fired, not that the compound is free of interference. Aggregation, redox cycling, metal
contamination, fluorescence and compound decomposition all produce false positives with no alert at
all, and only an assay can detect them.
