# Filter presets, and what each rule was actually derived from

Every rule below was derived from a specific dataset for a specific question. Applying one outside
that question is the most common way screening cascades go wrong. Cite the purpose, not just the
numbers.

## The published rule sets

| Rule | Criteria | Derived from | What it predicts |
|---|---|---|---|
| **Lipinski Ro5** (1997) | MW <= 500, cLogP <= 5, HBD <= 5, HBA <= 10; one violation tolerated | compounds that reached phase II | passive oral absorption - *not* activity, not safety |
| **Veber** (2002) | rotatable bonds <= 10, TPSA <= 140 | rat oral bioavailability, ~1100 compounds | oral bioavailability |
| **Egan** (2000) | TPSA <= 131.6, cLogP <= 5.88 | passive intestinal absorption | absorption |
| **Ghose** (1999) | MW 160-480, cLogP -0.4-5.6, 20-70 heavy atoms, MR 40-130 | CMC drug database | general drug-likeness; excludes small fragments by construction |
| **Rule of 3** (2003) | MW <= 300, cLogP <= 3, HBD <= 3, HBA <= 3, RotB <= 3 | fragment libraries | suitability as a *fragment starting point* |
| **Lead-likeness** | MW <= 350, cLogP <= 3.5, RotB <= 7 | lead-to-drug property drift | room to optimize without ending past Ro5 |

Two derived scores, both 0-1 or open-ended rather than pass/fail:

- **QED** (Bickerton 2012) - weighted desirability across eight properties. Useful for ranking within
  a set; meaningless as an absolute cut. Many approved drugs score below 0.4.
- **SA score** (Ertl & Schuffenhauer 2009) - 1 (easy) to 10 (hard), from fragment frequency in known
  synthesis. Above ~6 means "ask a chemist", not "impossible". Taxol scores ~5.9.

## Presets by purpose

**HTS hit triage (lead-like)**
```json
{"mw": [250, 400], "clogp": [1, 3.5], "rotatable_bonds": [null, 7], "tpsa": [40, 120]}
```
Leads gain roughly 50-100 Da and 1 log unit of cLogP during optimization. Start at the Ro5 ceiling
and you finish beyond it.

**Fragment screening**
```json
{"mw": [null, 300], "clogp": [null, 3], "hbd": [null, 3], "hba": [null, 3], "heavy_atoms": [8, 22]}
```
Ro5 is meaningless here. Judge fragments by ligand efficiency, not potency.

**Oral drug-like selection**
```json
{"mw": [null, 500], "clogp": [null, 5], "tpsa": [null, 140], "rotatable_bonds": [null, 10]}
```
with `require_rules: ["Lipinski", "Veber"]`.

**CNS penetrant**
```json
{"mw": [null, 400], "tpsa": [null, 90], "hbd": [null, 2], "clogp": [1.5, 4], "formal_charge": [0, 0]}
```
For anything serious here, the CNS MPO score (Wager 2010) is the better instrument; these are the
crude single-property cuts.

**Developability flags worth checking even when not filtering**

- **3/75 rule** (Hughes 2008): cLogP > 3 *and* TPSA < 75 correlates with in-vivo toxicity risk.
  Compute both and flag the quadrant rather than excluding.
- **Aromatic ring count > 4** correlates with poor developability (Ritchie & Macdonald 2009).
- **Fraction Csp3 < 0.2** flags flat, aromatic-heavy chemotypes with solubility risk.

## Where these rules do not apply

Say so explicitly rather than filtering anyway: natural products, macrocycles, peptides and
peptidomimetics, PROTACs and molecular glues, covalent inhibitors, antibiotics (notoriously
Ro5-breaking), prodrugs, and anything given by a non-oral route. A quarter of orally administered
approved drugs violate at least one Ro5 criterion.
