---
name: screen-library
description: >
  This skill should be used when the user asks to "screen this library", "filter these compounds",
  "apply Lipinski / Ro5 / lead-like filters", "clean up this SDF", "how many compounds pass", "triage
  this deck", "standardize these SMILES", or runs /screen-library. Covers the property and rule
  cascade, deduplication, salt stripping and the attrition table for compound files (CSV, SMI, SDF).
metadata:
  version: "0.1.0"
---

# Screening a compound library

Filtering a deck is easy; filtering it *defensibly* is the job. A cascade that discards 90% of a
library usually has a mis-set threshold in it, and the compounds it quietly removed may be the
interesting ones. Always show the attrition stage by stage and let the chemist argue with it.

## Procedure

1. **Look at the file before filtering it.** Call `screen_library` with no filters (or
   `standardize_library`) to get the parse failures, duplicate count and property distribution. Report:
   how many rows, how many parsed, how many duplicates, and the median and range of MW, cLogP, TPSA.
   A deck whose median MW is 230 is a fragment set and must not be judged by Ro5.
2. **Standardize first, and say so.** Salts, solvates and charge states wreck property calculations:
   a hydrochloride salt carries the counter-ion's mass into MW. `standardize_library` strips to the
   parent, neutralizes, and deduplicates by InChIKey. Mention how many compounds this changed.
3. **Choose thresholds from the purpose, not from habit.** Ask or infer which of these the user is
   doing, and set filters accordingly:
   - **Lead discovery / HTS triage** - lead-like: MW 250-400, cLogP 1-3.5, rotatable bonds <= 7.
     Leads grow during optimization, so starting at the Ro5 ceiling guarantees ending past it.
   - **Fragment screening** - Rule of 3: MW <= 300, cLogP <= 3, HBD/HBA <= 3. Never apply Ro5 here.
   - **Oral drug-like final selection** - Lipinski + Veber: MW <= 500, cLogP <= 5, HBD <= 5,
     HBA <= 10, rotatable bonds <= 10, TPSA <= 140.
   - **CNS target** - tighter: MW <= 400, TPSA <= 90, HBD <= 2, cLogP 2-4.
   - **Covalent or targeted-degrader work** - most property filters are inappropriate; say so rather
     than applying them.
   Read `references/filter-presets.md` for the numbers and what each rule was actually derived from.
4. **Run the cascade.** `screen_library` with `filters`, `require_rules`, and `alert_action`.
   Keep `alert_action="flag"` by default: see step 6.
5. **Report the attrition table, not just the survivors.** For each stage give in -> out and the
   most common failure reasons. If the pass rate is under ~5%, stop and question the thresholds
   before presenting the output; the tool flags this but the judgment is yours to state.
6. **Treat alerts as flags, never as an automatic cut.** Running this plugin's own example set with
   `alert_action="reject"` discards gefitinib and erlotinib - two marketed EGFR inhibitors - on
   Brenk's `Aliphatic_long_chain` and `triple_bond`. Report alerts in a column, let the chemist
   decide. Use `/alerts-review` when they want them interpreted.
7. **Write the output and hand back files.** `out_passed` and `out_failed` take .csv or .sdf; the
   failed file matters as much as the passed one, because that is what gets re-examined. Input
   columns (IDs, activity, vendor, catalogue numbers) are carried through.
8. **Offer the next step.** A filtered list is rarely the deliverable. Cluster it (`/hit-triage`) to
   see how many distinct chemotypes survived, and depict the top compounds so they can be looked at.

## Reporting

Give the chemist, in this order: the deck's shape before filtering, the cascade table, the pass rate,
the property distribution of what survived, the file paths, and the compounds that failed on a single
near-miss criterion (MW 505, cLogP 5.1) - those deserve a human look rather than silent loss.

Never describe a filtered set as "drug-like compounds" or imply the survivors are active. They are
compounds that satisfy numeric criteria. Activity is an assay result, not a descriptor.

## Common mistakes to avoid

- Applying Ro5 to a fragment deck, a natural-product set, a macrocycle series or a PROTAC set. All
  four break it by construction.
- Filtering before standardizing, so salt forms fail on MW.
- Deduplicating by SMILES string instead of InChIKey; two spellings of the same molecule are the same
  molecule, which is why `ASPIRIN` and `ASPIRIN_DUP` collapse in the example file.
- Reporting a pass count with no denominator or stage breakdown.
- Treating QED or synthetic accessibility as a verdict. SA score above ~6 says "ask a synthetic
  chemist", not "impossible"; taxol scores 5.9 and is made at scale.
