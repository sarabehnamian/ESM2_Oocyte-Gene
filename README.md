# Do protein language models work on the genes where IVF fails?

Real data throughout. ClinVar, gnomAD, UniProt, ESM-2. No simulation
anywhere in this pipeline.

## The question

Oocyte and IVF-failure genes are maternal-effect genes expressed only in
the egg and early embryo. TUBB8 is primate-specific. Their phenotype is
mostly invisible outside a fertility clinic: eggs that never mature, eggs
that never fertilise, embryos that arrest.

Protein language models learn from conservation across species. Genes
with little evolutionary depth may offer them little to learn from. So
the models could be weakest exactly where IVF fails.

The current fertility genetics literature mostly still annotates these
variants with SIFT and PolyPhen-2, from 2001 and 2010.

## Setup

```powershell
python -m pip install torch numpy matplotlib transformers biopython
```

`biopython` is only needed for step 08.

## Run order

```powershell
python 00_fetch_clinvar.py
python 01_fetch_gnomad.py
python 02_build_dataset.py
python 03_score_esm.py
python 04_evaluate.py
python 05_figures.py

# follow-up analyses
python 06_permutation_test.py      # runs from results/, no downloads
python 07_review_status.py         # runs from results/, no downloads
python 08_tubb8_paralog_positions.py   # needs rest.uniprot.org
python 09_pseudo_perplexity.py     # needs ESM-2 weights, GPU advised
python 10_alphamissense.py --am-file AlphaMissense_aa_substitutions.tsv.gz
python 11_mave_controls.py         # needs MaveDB and ESM-2 weights
python 12_figures_extra.py
```

Step 00 downloads about 440 MB once and caches it. Step 03 downloads the
ESM-2 weights, about 2.5 GB for the default 650M model.

Quick wiring check with a tiny model first:

```powershell
python 03_score_esm.py --model facebook/esm2_t6_8M_UR50D
```

## What each step does

| Step | Purpose | Writes |
|---|---|---|
| `00_fetch_clinvar.py` | Pathogenic missense variants, the positives | `data/clinvar_pathogenic.jsonl` |
| `01_fetch_gnomad.py` | Population variants above an allele frequency floor, the negatives | `data/gnomad_benign.jsonl` |
| `02_build_dataset.py` | UniProt sequences, position validation, merge | `data/variants.jsonl` |
| `03_score_esm.py` | ESM-2 zero-shot scores plus a BLOSUM62 baseline | `results/scores.jsonl` |
| `04_evaluate.py` | AUROC per gene and per panel, bootstrap intervals | `results/*.jsonl` |
| `05_figures.py` | Figures | `figures/*.png` |
| `06_permutation_test.py` | Exact permutation test of the panel gap over all C(14,4)=1001 label assignments, for ESM-2, BLOSUM62 and the ESM-minus-BLOSUM advantage; leave-one-out; min-AF 1e-3 | `results/permutation_*.jsonl` |
| `07_review_status.py` | ClinVar review status: 2+ star share per panel, AUROC by star class, panel test on 0-1 star positives only | `results/review_status_*.jsonl` |
| `08_tubb8_paralog_positions.py` | TUBB8 positions shared with the other human beta-tubulins vs TUBB8-specific; AUROC per class | `data/tubulin_paralogs.jsonl`, `results/tubb8_*.jsonl` |
| `09_pseudo_perplexity.py` | Per-gene ESM-2 pseudo-perplexity and marginal entropy, correlated with AUROC and advantage; random-init control | `results/pseudo_perplexity*.jsonl` |
| `10_alphamissense.py` | AlphaMissense as a third scorer | `data/alphamissense_panel.tsv.gz`, `results/*_am.jsonl` |
| `11_mave_controls.py` | Control-gene AUROC against deep mutational scanning labels from MaveDB | `results/mave_*.jsonl` |
| `12_figures_extra.py` | Figures 4-6 | `figures/fig4-6*.png` |

`common.py` holds the gene panels and metrics. Not a pipeline step.
Helpers for steps 06 onwards are appended at its end; nothing above them
changed. `manuscript_revision_notes.md` lists the text changes for the
docx that follow from steps 06 and 07.

## Why the negatives come from gnomAD

TUBB8 has 42 pathogenic missense variants in ClinVar and **zero benign
ones**. Several other panel genes are equally lopsided. There is no
AUROC without negatives.

So the negative set is missense variants carried by healthy people in
gnomAD above an allele frequency floor, default 1e-4. This is standard in
variant effect predictor benchmarking. It gives TUBB8 41 gnomAD
variants; 26 of them survive sequence validation and conflict removal,
against 38 positives.

## This pipeline has been run

Steps 00-05 were run end to end. The committed `results/` match the
650M numbers in the manuscript draft: ESM-2 gap 0.284 [0.153, 0.420],
BLOSUM62 gap 0.109 [-0.034, 0.250].

- 2,339 pathogenic variants from ClinVar
- 1,043 population variants from gnomAD
- 3,277 variants surviving sequence validation (2,321 pathogenic, 956 population)
- 14 of 28 genes with at least 5 variants in each class
- 4 oocyte genes, 10 control genes in the summary

`results/scores.jsonl` does not record whether wt-marginals or
masked-marginals produced it. Note this down when rerunning step 03.

## Results of the follow-up analyses run so far

Steps 06 and 07 run from the committed results with no downloads.

**06, exact permutation test** (`results/permutation_test.jsonl`).
Of the 1001 ways to place the 4 oocyte labels on the 14 included genes:

| quantity | oocyte | control | gap | assignments ≥ observed | p |
|---|---|---|---|---|---|
| ESM-2 AUROC | 0.619 | 0.903 | 0.284 | 1 / 1001 | 0.001 |
| BLOSUM62 AUROC | 0.598 | 0.707 | 0.109 | 47 / 1001 | 0.047 |
| ESM-2 minus BLOSUM62 | 0.021 | 0.196 | 0.175 | 2 / 1001 | 0.002 |

BLOSUM62 also shows some panel difference. The advantage gap is the
primary test. Leaving any one oocyte gene out keeps ESM and advantage
p ≤ 2/286 (`results/permutation_leave_one_out.jsonl`).

At min-AF 1e-3 (negatives filtered from the committed set): 3 oocyte vs
6 control genes, 84 assignments. ESM gap 0.153, p = 2/84; BLOSUM gap
0.027, p = 30/84; advantage gap 0.126, p = 6/84 = 0.071. This filter is
not identical to a rerun of steps 01-04 at 1e-3, because 13 positives
that were removed as conflicts at 1e-4 would come back; they are listed
in the output.

**07, review status** (`results/review_status_summary.jsonl`). 72.8% of
control positives have 2+ ClinVar stars, against 1.7% (1 of 58) of oocyte
positives. Using only 0-1 star positives on both sides, all 14 genes
still qualify: ESM gap 0.283 (p = 1/1001), BLOSUM gap 0.121 (35/1001),
advantage gap 0.162 (12/1001).

## Read before writing anything up

**Four oocyte genes is thin.** Only TUBB8 (38 vs 26), PATL2, PADI6 and
NLRP7 clear the threshold. Everything else has one to three pathogenic
variants. TUBB8 carries the argument almost by itself.

**gnomAD common does not mean harmless.** For recessive conditions a
pathogenic allele can sit at appreciable frequency. Try `--min-af 1e-3`
and check the conclusion survives.

**Ascertainment differs between panels.** Control gene variants come from
decades of clinical sequencing; oocyte gene variants come from small
infertility cohorts. That difference alone could produce a gap. It is the
strongest objection to this design and there is no clean fix within this
data.

**Sequence validation drops variants.** TUBB8 loses 13 percent, PATL2 25
percent. Every ClinVar positive validates; all 69 dropped variants are
gnomAD negatives, so the mismatch comes from the transcript behind
gnomAD's protein annotation, not from ClinVar. Step 02 reports the rate per gene. Anything above 50 percent
is flagged and should be excluded or isoform-mapped.

**If the gap crosses zero, the hypothesis fails.** Report that. A clean
negative on a well-posed question is worth more than a positive that
depended on picking the 8M model.

## Run locally: steps 08-11

The environment these steps were written in could not reach UniProt,
Hugging Face, Zenodo or MaveDB, so they have been tested only for code
paths (on synthetic inputs, nothing committed) and need running on your
machine. Each writes JSONL to `results/`; commit those files, then rerun
`python 12_figures_extra.py`.

```powershell
python -m pip install torch numpy matplotlib transformers biopython
```

**08, TUBB8 paralog positions.** Downloads 9 UniProt entries (a few KB).
Under 2 minutes on a CPU.

```powershell
python 08_tubb8_paralog_positions.py
```

Uses Biopython pairwise alignment to TUBB8 by default. If MAFFT is
installed, `--aligner mafft` uses it instead (this branch is untested).

**09, pseudo-perplexity.** ESM-2 650M weights, about 2.5 GB, the same
download as step 03 (cached by Hugging Face). One masked forward pass
per residue: 16,752 residues for the 14 included genes. Roughly 15-40
minutes per run on a recent GPU, most of a day on a CPU. Run it twice:
the trained model, then the random-init control.

```powershell
python 09_pseudo_perplexity.py --device cuda --batch-size 16
python 09_pseudo_perplexity.py --device cuda --batch-size 16 --random-init
```

Interrupted runs resume from `results/pseudo_perplexity.jsonl`. A
smoke test with a randomly initialised model already gave a sizeable
correlation between pseudo-perplexity and AUROC, so compare the trained
result with the random-init one before reading anything into it.

**10, AlphaMissense.** Download `AlphaMissense_aa_substitutions.tsv.gz`
from https://zenodo.org/records/8208688 (about 1.2 GB; check the size on
the record page). This is the protein-level table keyed by UniProt
accession, so it does not depend on the genome build; the hg19/hg38
files are not needed. Filtering takes about 5-10 minutes; the rest is
seconds.

```powershell
python 10_alphamissense.py --am-file C:\path\to\AlphaMissense_aa_substitutions.tsv.gz
```

This writes `data/alphamissense_panel.tsv.gz` (a few MB, the 28 panel
proteins only). Commit it; later runs need no argument.

**11, deep mutational scanning controls.** MaveDB downloads are a few
MB. ESM-2 scoring of the MAVE variants takes minutes on a CPU with
wt-marginals.

`data/mave_config.json` was written without access to MaveDB. Every
entry is marked `"verified": false`, and only BRCA1 (Findlay et al.
2018) has thresholds filled in. For each gene:

```powershell
python 11_mave_controls.py --search BRCA1          # list score sets
python 11_mave_controls.py --inspect urn:mavedb:...   # columns, score distribution
```

then set `urn`, `rule`, the thresholds from the paper, and
`"verified": true` in the config, and run

```powershell
python 11_mave_controls.py
```

The step refuses unverified entries (override with `--allow-unverified`)
and stops a gene whose wild-type match rate against UniProt is below 90%
(set `offset`).

## Novelty

Checked, twice, before this was written. The fertility genetics papers on
these genes use SIFT and PolyPhen-2, not modern protein language models,
and I found no evaluation of language models on this gene class. That is
two searches, not twenty. AlphaMissense scores every human missense
variant, so these genes may sit inside some large sweep. Worth one more
look before submission. Step 10 adds AlphaMissense as a scorer.
