# LTEE analysis (Supplementary Figure 4)

Loss-of-function (LoF) mutations in the 12 populations of the Long-Term Evolution Experiment (LTEE) compared with the
glucose fitness categories of the RB-TnSeq data (Price et al. 2018). Revised analysis for Reviewer 1 (major comment 6)
and Reviewer 2 (major comment 3).

## Run

Python 3 with pandas, numpy, scipy, matplotlib and openpyxl. From this folder:

```bash
python a09_ltee.py plus     # Ara+1 to Ara+6
python a09_ltee.py minus    # Ara-1 to Ara-6
python a09_ltee.py all      # the 12 populations
python a09b_deletion_audit.py
python a09c_match_barrick.py
python a10_fig_S4.py        # needs the three a09 runs
```

Each `a09_ltee.py` run takes a few minutes (10,000 permutations per subset). Results go to `Output/`.

| Script | What it does |
|---|---|
| `common.py` | Fitness loading (mean of the two replicates, `../Input_Data/fit_organism_Keio.tsv`) and category rule |
| `a09_ltee.py` | Parses the mutations, classifies LoF, removes amplification artefacts and runs the tests |
| `a09b_deletion_audit.py` | Audit of the large (≥ 5 kb) deletions of the Ara+ export |
| `a09c_match_barrick.py` | Matches every ALEdb clone to its curated GenomeDiff (population and mutator status) |
| `a10_fig_S4.py` | Supplementary Figure 4 (`Output/Fig_S4.svg/pdf/png`) |

## Data

| File | Source |
|---|---|
| `Input_Data/Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv` | ALEdb (Phaneuf et al. 2019), experiment "LTEE ARA" (Tenaillon et al. 2016): Ara+1 to Ara+6, 106 clones |
| `Input_Data/ALEdb_LTEE_AraMinus_mut.csv` | ALEdb, experiment "LTEE": Ara−1 to Ara−6, 96 clones |
| `Input_Data/fit_t_glucose_Price2018.tsv` | t scores of the two glucose experiments (set1IT003, set1IT004) from `fit_t.tab`, Price et al. 2018 supplementary data (https://genomics.lbl.gov/supplemental/bigfit/html/Keio/) |
| `Input_Data/enhancing_candidates_108_core.csv` | The 108 recurrent fitness-enhancing genes (fitness > 0.5 in ≥ 10 of 46 conditions); `significant_core` = also t > 4 in ≥ 1 of them (71 genes) |
| `LTEE-clone-curated/` | Curated breseq GenomeDiff files from https://github.com/barricklab/LTEE-Ecoli (commit 2307710, 23 Jun 2026; CC0, see `LICENSE` in that folder) |

## Methods in brief

- **LoF:** nonsense SNPs; coding indels whose length is not a multiple of 3; IS insertions inside a coding region;
  deletions that remove whole or partial genes. Missense, synonymous, in-frame and intergenic mutations are not LoF.
- **Amplification artefacts:** in the curated GenomeDiffs, eight deletion rows carry `within=<AMP>:<copy>`, i.e. they
  lie inside one copy of an amplification and no gene is lost. ALEdb drops this attribute, so these rows are removed
  (`a09b`, `a09c`). They include the only apparent LoF of *rpoS*.
- **Universe:** genes with glucose fitness that are in the EcoCyc K-12 gene table. REL606 gene names are matched to
  K-12 b-numbers by name; IS elements and `ECB_` loci are not mapped.
- **Mutators:** ALEdb only includes clones without elevated point-mutation rates. IS-mutator clones (Ara+1 from
  5,000 generations; one Ara−5 clone at 30,000) are taken from the GenomeDiff metadata and analysed separately.
- **Tests:** hypergeometric test per fitness category; event-level permutation (10,000 random gene sets, single-gene
  events drawn in proportion to gene length, multi-gene deletions as blocks of the same number of consecutive genes);
  Mann-Whitney test of glucose fitness. Subsets: all LoF, point-like LoF, parallel (≥ 2 populations), persistent
  (≥ 2 time points) and without IS-mutator clones.
