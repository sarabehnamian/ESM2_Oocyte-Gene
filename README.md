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
python -m pip install torch numpy matplotlib transformers
```

## Run order

```powershell
python 00_fetch_clinvar.py
python 01_fetch_gnomad.py
python 02_build_dataset.py
python 03_score_esm.py
python 04_evaluate.py
python 05_figures.py
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

`common.py` holds the gene panels and metrics. Not a pipeline step.

## Why the negatives come from gnomAD

TUBB8 has 42 pathogenic missense variants in ClinVar and **zero benign
ones**. Several other panel genes are equally lopsided. There is no
AUROC without negatives.

So the negative set is missense variants carried by healthy people in
gnomAD above an allele frequency floor, default 1e-4. This is standard in
variant effect predictor benchmarking. It gives TUBB8 41 negatives
against its positives.

## This pipeline has been run

Every step was executed end to end before it was handed over, using
**ESM-2 8M**, the smallest model, purely to check the wiring. Real
numbers from that run:

- 2,339 pathogenic variants from ClinVar
- 1,043 population variants from gnomAD
- 3,277 variants surviving sequence validation
- 14 of 28 genes with at least 5 variants in each class
- 4 oocyte genes, 10 control genes in the summary

The 8M smoke test gave an ESM gap of 0.202 with interval [0.026, 0.430],
and BLOSUM62 crossing zero at 0.109 [-0.027, 0.242].

**Do not report those numbers.** An 8M model is not a result. Rerun with
the default 650M model, which is what step 03 uses if you pass nothing.
The direction may change.

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
percent, from UniProt canonical disagreeing with the RefSeq transcript
ClinVar used. Step 02 reports the rate per gene. Anything above 50 percent
is flagged and should be excluded or isoform-mapped.

**If the gap crosses zero, the hypothesis fails.** Report that. A clean
negative on a well-posed question is worth more than a positive that
depended on picking the 8M model.

## Novelty

Checked, twice, before this was written. The fertility genetics papers on
these genes use SIFT and PolyPhen-2, not modern protein language models,
and I found no evaluation of language models on this gene class. That is
two searches, not twenty. AlphaMissense scores every human missense
variant, so these genes may sit inside some large sweep. Worth one more
look before submission.
