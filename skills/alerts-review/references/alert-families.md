# What the common alert families actually mean

Grouped by mechanism, because the mechanism decides the experiment. Names in `code` font are the
descriptions RDKit's catalogs return.

## Assay interference (the PAINS core)

| Family | Mechanism | Settling experiment |
|---|---|---|
| Rhodanines, ene-rhodanines (`ene_rhod_A`) | Michael acceptor, photoreactivity, aggregation; hit in a very high fraction of screens | orthogonal non-fluorescent readout; SAR across close analogues |
| Catechols (`catechol_A`) | redox cycling, oxidation to quinone, metal chelation | assay with and without reducing agent; check for H2O2 production |
| Quinones | redox cycling, ROS generation, covalent addition | catalase control; thiol-reactivity test |
| Phenolic Mannich bases (`mannich_A`) | eliminate to a quinone methide electrophile | mass spectrometry for protein adducts |
| Acylhydrazones, hydroxyphenylhydrazones | hydrolysis to reactive aldehyde, metal chelation, photochemistry | fresh vs aged DMSO stock; metal-free buffer |
| beta-Diketones / curcuminoids | metal chelation, Michael addition, membrane perturbation, fluorescence | orthogonal format; counter-screen unrelated targets |
| 2-Aminothiazoles, aminothiophenes | frequent hitters across unrelated targets | hit rate across the user's own historical screens |

## Reactivity (electrophiles)

Acyl halides, sulfonyl halides, epoxides, aziridines, alkyl halides, aldehydes, maleimides,
acrylamides, vinyl sulfones, isocyanates, beta-lactams.

These react with protein nucleophiles, with buffer components, or with the compound's own storage
solvent. **Context decides everything**: an acrylamide in a covalent-inhibitor programme is the
designed warhead, and beta-lactams are an antibiotic class, not a liability. Flag reactivity, ask
whether it is intended, and only then advise. For unintended reactivity, test stability in buffer and
in glutathione.

## Metabolic and toxicological liabilities

| Family | Concern |
|---|---|
| Nitroaromatics | nitroreduction to reactive intermediates; frequent Ames positives |
| Anilines | CYP oxidation to hydroxylamines and quinone imines; idiosyncratic tox |
| Thiophenes, furans | CYP bioactivation to S-oxides and epoxides |
| Hydrazines, azo groups | genotoxicity, metabolic instability |
| Long aliphatic chains (`Aliphatic_long_chain`) | promiscuity, aggregation, poor solubility and properties |
| Polyhalogenated aromatics | persistence, lipophilicity, tox |

These are *development* risks. They belong in a discussion about a lead series, not in a filter that
runs before anyone has looked at the data. Gefitinib carries `Aliphatic_long_chain`; erlotinib
carries `triple_bond`; both are marketed drugs.

## What alerts cannot detect

No catalog flags these, and each produces false positives with a perfectly clean structure:

- **Colloidal aggregation** - the single largest source of HTS false positives. Detect with 0.01%
  Triton X-100 or with a detergent-sensitivity counter-screen, and by checking for a steep Hill
  slope and enzyme-concentration dependence.
- **Compound decomposition or impurity** - the activity belongs to something else in the well. Test
  fresh solid and check purity by LC-MS.
- **Fluorescence or absorbance interference** - the compound is read as signal. Use an orthogonal
  readout format.
- **Metal contamination** from synthesis, giving metal-dependent inhibition.
- **Redox cycling** producing H2O2 in the assay buffer; test with catalase.

A compound with no alerts is a compound with no *known* problem substructure. That is a weaker
statement than "clean", and worth saying out loud.

## Sources

Baell & Holloway, *J Med Chem* 53 (2010) 2719 - the PAINS filters. Brenk et al., *ChemMedChem* 3
(2008) 435 - unwanted functionality. Baell & Nissink, *ACS Chem Biol* 13 (2018) 36 - a correction to
the over-application of PAINS filters, worth reading before using them as a blacklist.
