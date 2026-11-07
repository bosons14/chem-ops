---
name: molecule-profile
description: >
  This skill should be used when the user asks to "profile this compound", "what are the properties
  of this SMILES", "is this drug-like", "calculate logP/TPSA/MW", "compute ligand efficiency",
  "what's the LE/LLE of this hit", "compare these two compounds", or runs /molecule-profile. Covers
  descriptors, drug-likeness rule sets, synthetic accessibility and efficiency metrics for one
  compound or a small set.
metadata:
  version: "0.1.0"
---

# Profiling a compound

A chemist asking "is this drug-like?" wants a reading of where the compound sits and what would need
fixing, not a pass/fail stamp. Give the numbers, name the ones that are out of range, and say what
they imply for the next design cycle.

## Procedure

1. **Profile it.** `molecule_profile` with the SMILES. For a single compound set
   `canonical_tautomer=true` - it is slow in bulk but correct, and tautomer choice changes HBD/HBA
   counts. Pass `image_path` so the structure is actually shown with alert atoms highlighted.
2. **Add potency if the user has any.** `potency` with `potency_kind` (`pic50`/`pki`/`pkd`, or
   `ic50_nm`/`ki_nm`/`kd_nm`) gives ligand efficiency metrics. These are the most useful numbers in
   the whole profile during optimization:
   - **LE** (kcal/mol per heavy atom) >= 0.3 is the usual target. A large weak compound is a worse
     starting point than a small moderate one.
   - **LLE** (pActivity - cLogP) >= 5 is the usual target. LLE is what separates real affinity from
     potency bought with grease.
   - **LELP** (cLogP / LE) above ~10 says the potency came from lipophilicity.
   Imatinib, for scale, has LE ~0.29 and LELP ~16 - a reminder that approved drugs are not always
   efficient, and that these are optimization guides rather than gates.
3. **Interpret the descriptors in context.** Report what matters for the user's stage:
   - MW and cLogP against the project's room to grow, not against 500 and 5 in the abstract.
   - TPSA for permeability (< 90 for CNS, < 140 for oral absorption).
   - HBD count - the single property most predictive of poor permeability; each one costs.
   - Rotatable bonds for conformational entropy on binding and for oral bioavailability.
   - Fraction Csp3 and aromatic ring count for developability and solubility risk.
   - Unspecified stereocentres - flag these loudly. A compound drawn without stereochemistry that
     has three centres is eight compounds, and the assay result may belong to one of them.
4. **Name the failed criteria specifically.** "Fails Lipinski" is useless; "MW 853 and 14 HBA, so
   two Ro5 violations" tells the chemist which handle to pull.
5. **Comparing compounds?** Profile them together and present a table of the deltas, with the
   property that differs most called out. If they share a core, use `series_core` and depict them
   aligned so the comparison is visual as well as numeric.

## What to say and what not to say

Say: where the compound sits, which properties are out of the usual range for this purpose, what the
efficiency metrics suggest about the series, and what a chemist might change.

Do not say: that a compound is active, selective, safe, permeable or orally bioavailable. None of
those are computed here. cLogP is a *calculated* estimate that routinely differs from measured logD
by a log unit or more; call it cLogP, not logP, and never present it as measured. Solubility,
metabolic stability, hERG and tox are not computed by this plugin at all - if the user asks, say
they need measured data or a dedicated predictive model.

If the compound profile feeds a decision about a real programme, say plainly that this is a
computational first pass for a chemist to judge.
