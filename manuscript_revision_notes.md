# Manuscript revision notes

Changes for `ESM2_Oocyte-Gene_Benchmark_Methods-Results.docx`. The docx
itself has not been edited. Each item gives the current text (**Before**)
and replacement text (**After**) to paste. "Before" snippets are copied
from the docx; where only part of a paragraph changes, only that part is
quoted.

Every number in the "After" text comes from a committed file, named in
brackets after the item. Items marked **[needs local run]** wait on
steps 08-11. Leave them out or keep them as placeholders until those
steps have run.

---

## 0. Numbers that do not trace to a committed file

Fix these before anything else. Each one is either wrong or comes from a
run whose output was never saved.

| Docx statement | Committed files say | Action |
|---|---|---|
| Leave-one-out gaps 0.267 (NLRP7), 0.280 (PATL2), 0.293 (PADI6), 0.302 (TUBB8) | 0.270, 0.281, 0.290, 0.294 (`results/permutation_leave_one_out.jsonl`, from `results/auroc_by_gene.jsonl`) | Replace, see 4.6 |
| wt-marginals gap 0.286 [0.149, 0.435] | No committed file. `results/scores.jsonl` does not record which method produced it | Rerun step 03 with each method into separate copies and commit both, or drop the sentence |
| "masked-marginals ... primary analysis" | Step 03 defaults to wt-marginals, and the README run order uses the default | Confirm which method produced `results/scores.jsonl` and say so in Methods |
| min-AF 1e-3: negatives fall "from 1,043 to 299" | `data/gnomad_benign.jsonl` has 298 variants with AF ≥ 1e-3. Filtering the scoring set gives 268 negatives after validation and conflict removal | See 4.6 |
| min-AF 1e-3: gap 0.160 [0.004, 0.328]; advantage +0.152 vs +0.030 | Filtering the committed scoring set gives an ESM gap of 0.153 and advantages of +0.163 vs +0.037 (`results/permutation_test.jsonl`, setting `min_af=0.001`). The difference: a true rerun adds back 13 positives that were removed as conflicts at 1e-4 and never scored | Rerun steps 01-04 at 1e-3 and commit the outputs, or use the filtered numbers with the caveat in 4.6 |
| TUBB8 "42 pathogenic ... zero benign" (Methods 00) and "41 for TUBB8" (Methods 01), versus 38/26 in Table 1 | Both are correct but measure different stages: 42 unique ClinVar positives → 38 after removing 4 that also appear in gnomAD. 41 gnomAD variants → 30 after sequence validation (11 wild-type mismatches) → 26 after the same 4 conflicts (`data/clinvar_pathogenic.jsonl`, `data/gnomad_benign.jsonl`, `data/counts_dataset.jsonl`) | State the path once, see 2.2 |
| Drop rates "13 percent for TUBB8 and 25 percent for PATL2", explained as UniProt vs the RefSeq transcript ClinVar used | All 2,339 ClinVar positives validate. All 69 dropped variants are gnomAD negatives (TUBB8 11 of 41, PATL2 7 of 22, PTEN 13 of 19, TP53 10 of 21, LDLR 13 of 39) | Correct the explanation, see 2.3 |
| "all four oocyte genes ... with intervals that either touch or cross the chance line" | TUBB8's interval is [0.509, 0.775], which does not touch 0.5 | See 4.2 |
| "a six-step pipeline of numbered Python scripts sharing two modules" | There is one shared module, `common.py` | See 2.1 |

---

## 1. Title and abstract

The title can stay. If an abstract is added, lead with the advantage-gap
permutation test (section 4.3), not the bootstrap interval.

---

## 2. Methods

### 2.1 Pipeline overview

**Before**
> The analysis is implemented as a six-step pipeline of numbered Python scripts sharing two modules.

**After**
> The analysis is implemented as numbered Python scripts sharing one module of helpers. Steps 00 to 05 build the dataset, score it and produce the main figures; steps 06 to 12 add the permutation test, the review-status check, the TUBB8 paralog analysis, pseudo-perplexity, AlphaMissense, deep mutational scanning controls and the extra figures.

### 2.2 00_fetch_clinvar.py, the TUBB8 count

**Before**
> TUBB8, the flagship oocyte maturation arrest gene, has 42 pathogenic missense variants in ClinVar and zero benign ones.

**After**
> TUBB8, the flagship oocyte maturation arrest gene, has 42 distinct pathogenic missense variants in ClinVar and no benign ones. After the steps below, 38 of them remain in the scoring set; the other four are also present in gnomAD above the frequency floor and are removed as ambiguous.

[`data/clinvar_pathogenic.jsonl`, `data/counts_dataset.jsonl`]

**Before** (01_fetch_gnomad.py)
> This yields 1,043 population variants across the 28 genes, including 41 for TUBB8, which converts a gene with no usable negatives into one with a nearly balanced comparison.

**After**
> This yields 1,043 population variants across the 28 genes, including 41 for TUBB8. Of those 41, 11 fail sequence validation and 4 are also ClinVar pathogenic, leaving 26 TUBB8 negatives against 38 positives in the scoring set.

[`data/gnomad_benign.jsonl`, `data/counts_dataset.jsonl`]

### 2.3 02_build_dataset.py, which variants are dropped

**Before**
> Observed drop rates were low for most genes (zero for the majority) but reached 13 percent for TUBB8 and 25 percent for PATL2.

**After**
> Every ClinVar positive passed validation. All 69 dropped variants were gnomAD negatives, which points to the transcript used for the gnomAD protein annotation rather than to UniProt disagreeing with ClinVar. Drop rates were zero for 17 of the 28 genes and reached 13 percent for TUBB8 (11 of 41 gnomAD variants), 25 percent for PATL2 overall (7 of 22 gnomAD variants) and, among controls, 13 of 19 gnomAD variants for PTEN and 10 of 21 for TP53. Because validation only checks the wild-type residue, a gnomAD variant annotated on a shifted transcript can occasionally pass by chance; this would add noise to the negative class.

[`data/counts_dataset.jsonl`, `data/counts_gnomad.jsonl`]

### 2.4 03_score_esm.py, scoring method

**Before**
> Both were run and are compared in the Results.

**After** (fill in once confirmed)
> Both were run. The results reported here use [wt-marginals / masked-marginals]; the other method is compared in the Robustness section.

### 2.5 04_evaluate.py, the panel test

**Before**
> The difference between panels is estimated the same way, producing a direct interval on the control-minus-oocyte gap for each scorer, which is the quantity of interest.

**After**
> The difference between panels is estimated the same way, giving a descriptive interval on the control-minus-oocyte gap for each scorer. With four oocyte genes the coverage of a gene-level bootstrap interval is uncertain, so the panel difference is tested with an exact permutation test (06_permutation_test.py). Per-gene AUROCs are held fixed and the oocyte label is placed on every possible set of four of the 14 included genes, C(14,4) = 1,001 assignments. For each assignment the gap is the control mean minus the oocyte mean. The one-sided p-value is the fraction of assignments whose gap is at least the observed gap; the observed assignment is one of them, so the smallest attainable p is 1/1,001. The test is applied to three quantities: the ESM-2 gap, the BLOSUM62 gap, and the gap in ESM-2 advantage, where the advantage for a gene is its ESM-2 AUROC minus its BLOSUM62 AUROC. The advantage gap is the primary test, because it asks whether the language model's gain over a fixed substitution matrix is smaller on oocyte genes, which panel-wide differences in label quality would not by themselves produce.

### 2.6 New Methods paragraphs

Add after 05_figures.py:

> **06_permutation_test.py.** Exact permutation test as described above, repeated with each oocyte gene left out in turn (C(13,3) = 286 assignments) and with the stricter allele frequency floor of 1e-3.
>
> **07_review_status.py.** ClinVar review status was converted to gold stars (0 for no assertion criteria, 1 for single submitter, 2 for multiple submitters without conflict, 3 for expert panel, 4 for practice guideline). Per-gene AUROC was recomputed with positives split into 0-1 star and 2+ star, each against all of that gene's negatives, and the panel comparison was repeated using only 0-1 star positives in both panels, with the same inclusion rule (at least five variants per class).
>
> **08_tubb8_paralog_positions.py** [needs local run]. The nine human beta-tubulins (TUBB, TUBB1, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6, TUBB8) were taken from UniProt and each paralog was aligned to TUBB8 by global pairwise alignment (BLOSUM62, gap open -10, extend -0.5). A TUBB8 position was labelled shared if its residue was identical in all eight other paralogs, and TUBB8-specific otherwise; a threshold of six of eight was used as a sensitivity check. ESM-2 and BLOSUM62 AUROC were computed within each class, with class-stratified bootstrap intervals, and the difference between classes was tested by permuting class labels within pathogenic and within benign variants (10,000 permutations).
>
> **09_pseudo_perplexity.py** [needs local run]. For each included gene, every residue was masked in turn and the model's probability of the true residue recorded; pseudo-perplexity is the exponential of the negative mean log-probability. Proteins longer than 1,022 residues were scored in windows of 1,022 residues centred on each position, so for these genes the model never sees the whole protein. The mean entropy of the model's predicted distribution over the 20 standard amino acids, from one unmasked pass, was also recorded. Both were correlated with per-gene ESM-2 AUROC and with ESM-2 advantage (Spearman, two-sided permutation p, 100,000 permutations). The same computation with a randomly initialised model of the same architecture served as a control for effects of sequence composition.
>
> **10_alphamissense.py** [needs local run]. Precomputed AlphaMissense scores (Cheng et al. 2023) were matched to the scoring set by UniProt accession, wild-type residue, position and mutant residue. Per-gene AUROC for all three scorers was computed on the variants AlphaMissense covers, and the panel tests were repeated with AlphaMissense included.
>
> **11_mave_controls.py** [needs local run]. For BRCA1, TP53, PTEN and MSH2, one published deep mutational scanning score set per gene was taken from MaveDB and variants were labelled functional or non-functional using the thresholds of the original study [list each: score set URN, assay, threshold, citation]. Intermediate scores were excluded. ESM-2 scores were computed with the same log-likelihood ratio as step 03.

---

## 3. Figure captions

### 3.1 Figure 1

**Before**
> Figure 1. Per-gene ESM-2 AUROC with 95 percent bootstrap intervals. Counts are shown as pathogenic/benign. Oocyte genes in green, conserved controls in purple.

**After**
> Figure 1. Per-gene ESM-2 AUROC with 95 percent bootstrap intervals. Counts are shown as pathogenic/benign. Oocyte genes in green, conserved controls in purple. The dashed line marks chance (0.5).

(The "chance" label overlapping the title is fixed in `05_figures.py`; it now sits at the foot of the dashed line. Re-insert `figures/fig1_auroc_by_gene.png`. It was regenerated without Georgia installed, so rerun `python 05_figures.py` on a machine with Georgia if the font must match.)

### 3.2 Figure 2

**Before**
> Figure 2. Left: mean AUROC per panel for both scorers. Right: the control minus oocyte gap with 95 percent bootstrap intervals. ESM-2 loses 0.284; BLOSUM62 shows no detectable difference.

**After**
> Figure 2. Left: mean AUROC per panel for both scorers. Right: the control-minus-oocyte gap with 95 percent bootstrap intervals, shown as a description of uncertainty; the panel test is the exact permutation test in Figure 4. The ESM-2 gap is 0.284 and the BLOSUM62 gap is 0.109.

### 3.3 New Figure 4

> Figure 4. Exact permutation null distributions over all 1,001 ways of placing the four oocyte labels on the 14 included genes. Left: gap in ESM-2 advantage over BLOSUM62 (observed 0.175, p = 2/1,001). Centre: ESM-2 gap (observed 0.284, p = 1/1,001). Right: BLOSUM62 gap (observed 0.109, p = 47/1,001). Dark bars are assignments at least as extreme as observed.

[`figures/fig4_permutation_null.png`, `results/permutation_test.jsonl`]

### 3.4 Figures 5 and 6 [needs local run]

> Figure 5. TUBB8 ESM-2 and BLOSUM62 AUROC at positions shared with all eight other human beta-tubulins and at TUBB8-specific positions, with 95 percent bootstrap intervals.
>
> Figure 6. Per-gene ESM-2 AUROC (left) and ESM-2 advantage over BLOSUM62 (right) against ESM-2 pseudo-perplexity of the wild-type protein.

---

## 4. Results

### 4.1 Dataset composition

No number changes. Optionally add after the first paragraph:

> The panels also differ in how the positive labels were curated: 72.8 percent of control-gene positives (1,629 of 2,238) carry a ClinVar review status of two stars or more, against 1.7 percent of oocyte-gene positives (1 of 58). The consequences are examined below.

[`results/review_status_summary.jsonl`]

### 4.2 Per-gene discrimination

**Before**
> all ten control genes have ESM-2 AUROC above 0.76 and their intervals lie entirely above chance, whereas all four oocyte genes fall between 0.577 and 0.650 with intervals that either touch or cross the chance line.

**After**
> all ten control genes have ESM-2 AUROC above 0.76 and their intervals lie entirely above chance, whereas all four oocyte genes fall between 0.577 and 0.650. The intervals for PADI6, PATL2 and NLRP7 cross chance; the interval for TUBB8 lies just above it.

[`results/auroc_by_gene.jsonl`]

### 4.3 The panel gap: new headline test

**Before** (second paragraph of "The panel gap and why the baseline matters")
> The BLOSUM62 baseline discriminates between these explanations. It reaches 0.707 (0.639 to 0.773) on controls and 0.598 (0.476 to 0.733) on oocyte genes, a gap of 0.109 with an interval of [-0.034, 0.250] that crosses zero. A scorer carrying no gene-specific information at all shows no detectable panel difference.

**After**
> The BLOSUM62 baseline helps separate these explanations. It reaches 0.707 (0.639 to 0.773) on controls and 0.598 (0.476 to 0.733) on oocyte genes, a gap of 0.109. BLOSUM62 carries no gene-specific information, yet it also shows some panel difference: in the exact permutation test, 47 of the 1,001 possible assignments of the oocyte label give a BLOSUM62 gap at least this large (one-sided p = 0.047). Part of the ESM-2 gap is therefore shared with a scorer that knows nothing about the genes, and is most simply explained by differences between the panels in their variants or labels.
>
> The question is whether ESM-2 loses more than that. For each gene we took the ESM-2 advantage, its ESM-2 AUROC minus its BLOSUM62 AUROC. The mean advantage is 0.196 on control genes and 0.021 on oocyte genes, a gap of 0.175. Only 2 of the 1,001 label assignments give an advantage gap this large (one-sided p = 0.002). For comparison, the ESM-2 gap itself is 0.284, and no other assignment reaches it (p = 1/1,001, the smallest value the test can return). The language model's gain over the substitution matrix is much smaller on the oocyte genes, and this is the main result.

[`results/permutation_test.jsonl`, setting `min_af=0.0001`]

**Before** (first paragraph of the same section)
> The gap is 0.284 with an interval of [0.153, 0.420], which excludes zero.

**After**
> The gap is 0.284 (bootstrap interval 0.153 to 0.420, reported as a description of uncertainty; the test of the gap follows below).

### 4.4 The advantage paragraph

**Before**
> The clearest way to state the result is as the advantage that the language model provides over the substitution matrix. On control genes that advantage is +0.196 AUROC. On oocyte genes it is +0.021. A model with 650 million parameters, trained on hundreds of millions of protein sequences, delivers approximately two hundredths of an AUROC over a fixed substitution matrix from 1992 on the genes where IVF fails. Whatever ESM-2 has learned that makes it useful on BRCA1 and LDLR does not transfer to TUBB8 and NLRP7.

**After**
> On the four oocyte genes, ESM-2 (650 million parameters) gains 0.021 AUROC over BLOSUM62, a fixed substitution matrix from 1992. On the control genes it gains 0.196. Whatever ESM-2 adds over BLOSUM62 on BRCA1 and LDLR, it adds little of it on TUBB8, PADI6, PATL2 and NLRP7. On PADI6 its AUROC is below BLOSUM62's (0.637 against 0.748).

[`results/auroc_by_gene.jsonl`]

### 4.5 Mechanism paragraph (soften until 09 and 08 have run)

**Before**
> This is consistent with the mechanism proposed at the outset. Protein language models learn from patterns of conservation across evolutionary time. TUBB8 is primate-specific and expressed only in oocytes and early embryos; the oocyte panel more broadly consists of rapidly evolving, lineage-restricted maternal-effect genes. The evolutionary signal these models depend on is thin or absent for precisely this gene class, and the results are what that account predicts.

**After**
> One possible explanation is the one proposed at the outset: protein language models learn from conservation across evolutionary time, and many oocyte and maternal-effect genes evolve quickly or are restricted to particular lineages, so the signal these models rely on may be thin. These results do not test that explanation. TUBB8 is also a poor fit for it as stated. TUBB8 arose in primates, but it is a beta-tubulin, one of the most conserved protein families, and most of its residues are shared with the other human beta-tubulins that ESM-2 has seen in training. If thin evolutionary signal were the cause, ESM-2 should find TUBB8 hard to predict and should fail mainly at TUBB8-specific positions. Both predictions can be checked directly, and are examined in [section X, once steps 08 and 09 are run]. Other explanations remain open, including differences in how the oocyte-gene variants were ascertained and labelled, and pathogenic mechanisms in these genes (for example, dominant-negative effects in TUBB8) that a single-site likelihood ratio may not capture.

**After 09 has run**, replace the "These results do not test that explanation" sentence onward with the result, whichever way it goes:
> Across the 14 included genes, ESM-2 pseudo-perplexity correlated with per-gene ESM-2 AUROC at Spearman rho = [ ] (permutation p = [ ]) and with ESM-2 advantage at rho = [ ] (p = [ ]). A randomly initialised model of the same architecture gave rho = [ ] and [ ]. TUBB8 ranked [ ] of [ ] genes by pseudo-perplexity.

[`results/pseudo_perplexity_correlations.jsonl`]. If the trained-model correlation is no stronger than the random-init one, say the explanation is not supported by this test.

**After 08 has run**:
> Within TUBB8, ESM-2 AUROC was [ ] at the [ ] positions shared with all eight other beta-tubulins ([ ] pathogenic, [ ] benign) and [ ] at the [ ] TUBB8-specific positions ([ ], [ ]); the difference was [ ] (permutation p = [ ]).

[`results/tubb8_position_classes.jsonl`]

### 4.6 Robustness

**Before**
> First, leave-one-out across the four oocyte genes gives gaps of 0.267 (dropping NLRP7), 0.280 (PATL2), 0.293 (PADI6) and 0.302 (TUBB8). No single gene drives the result, and notably the gap widens rather than collapses when TUBB8, the best-powered gene, is removed.

**After**
> First, leave-one-out across the four oocyte genes gives ESM-2 gaps of 0.270 (dropping NLRP7), 0.281 (PATL2), 0.290 (PADI6) and 0.294 (TUBB8), and advantage gaps of 0.183, 0.181, 0.131 and 0.204. In each case the exact permutation test over the 286 remaining assignments gives p ≤ 2/286 for both the ESM-2 gap and the advantage gap. No single gene drives the result, and the gap does not shrink when TUBB8, the best-powered gene, is removed.

[`results/permutation_leave_one_out.jsonl`]

**Before**
> Third, the definition of the benign class was tightened from an allele frequency floor of 1e-4 to 1e-3. This preserves the direction, with a gap of 0.160, but the interval widens to [0.004, 0.328] and its lower bound nearly touches zero. The cause is mechanical: the stricter floor reduces the negative set from 1,043 to 299 variants, PTEN loses all its negatives, and MSH2, TP53 and RB1 fall below the inclusion threshold, leaving only three oocyte and six control genes. The advantage-over-baseline pattern is nonetheless preserved (+0.152 on controls versus +0.030 on oocyte genes). This check is reported as consistent with the primary analysis rather than as independent confirmation of it.

**After** (numbers from filtering the committed scoring set; replace with the outputs of a full rerun if you do one)
> Third, the benign class was restricted to gnomAD variants with allele frequency of at least 1e-3 instead of 1e-4. This leaves 268 negatives in the scoring set. PTEN loses all its negatives, and MSH2, TP53, RB1 and PATL2 fall below the inclusion threshold, leaving three oocyte genes (TUBB8, PADI6, NLRP7) and six control genes, so there are only 84 possible label assignments and the smallest attainable p is 1/84. The ESM-2 gap is 0.153 (p = 2/84 = 0.024). The BLOSUM62 gap falls to 0.027 (p = 30/84 = 0.36). The advantage gap is 0.126 (+0.163 on controls against +0.037 on oocyte genes; p = 6/84 = 0.071). The direction holds, but with this few genes the advantage gap is not significant on its own, and this check is consistent with the primary analysis rather than independent confirmation of it. These values come from filtering the 1e-4 scoring set; a full rerun at 1e-3 would also restore 13 ClinVar positives (8 in CFTR, 3 in NLRP7, 1 in LDLR, 1 in TUBB8) that were removed as conflicts at 1e-4.

[`results/permutation_test.jsonl`, setting `min_af=0.001`; the 13 variants are listed in its `unscored_positives` field]

### 4.7 New subsection: review status

Add after Robustness (or as a fourth robustness check):

> **ClinVar review status.** The positives in the two panels differ in how thoroughly they were reviewed: 72.8 percent of control-gene positives have two or more ClinVar stars, against one of 58 oocyte-gene positives (a PATL2 variant). If better-reviewed positives are easier to separate, the panel gap could reflect curation rather than biology. Two checks argue against this. Within the ten control genes, ESM-2 AUROC computed with 0-1 star positives only (mean 0.882) is close to that with 2+ star positives only (mean 0.902), each against the same negatives. And when the panel comparison is repeated with only 0-1 star positives on both sides (all 14 genes still qualify; PATL2 drops to 5 positives), the ESM-2 gap is 0.283 (p = 1/1,001), the BLOSUM62 gap is 0.121 (p = 35/1,001) and the advantage gap is 0.162 (p = 12/1,001). The result does not depend on the well-reviewed control variants.

[`results/review_status_summary.jsonl`, `results/review_status_by_gene.jsonl`]

### 4.8 New subsections, AlphaMissense and deep mutational scans [needs local run]

> **AlphaMissense.** [Match rate per gene. Per-gene AUROC. Panel means. Permutation p for the AlphaMissense gap and the AlphaMissense-minus-BLOSUM62 advantage gap.] Note in the text that AlphaMissense was trained with common human and primate variants as its benign class, so some of the gnomAD negatives here were likely benign-labelled during its training.
>
> [`results/auroc_by_gene_am.jsonl`, `results/panel_summary_am.jsonl`, `results/permutation_test_am.jsonl`, `results/alphamissense_match.jsonl`]
>
> **Labels independent of ClinVar.** [For each of BRCA1, TP53, PTEN and MSH2: score set, threshold, number of functional and non-functional variants, ESM-2 and BLOSUM62 AUROC against the assay labels, and the ClinVar-based AUROC for comparison.]
>
> [`results/mave_auroc.jsonl`]

---

## 5. Limitations

### 5.1 Ascertainment paragraph

**Before**
> The BLOSUM62 control addresses the most direct form of this objection, since a difference in label quality alone should degrade both scorers, but it does not exclude subtler forms, such as pathogenic variants in the oocyte panel being enriched for effects that sequence context cannot capture.

**After**
> The BLOSUM62 control addresses part of this objection, since a difference in label quality alone should degrade both scorers. BLOSUM62 does show a smaller panel difference of its own (p = 0.047), which is consistent with some such difference being present. The advantage-gap test asks whether ESM-2 loses more than that, but it does not exclude subtler explanations, such as pathogenic variants in the oocyte panel being enriched for effects that sequence context cannot capture.

### 5.2 New paragraph: circularity between ClinVar labels and the scorers

> The ClinVar labels are not independent of the kind of evidence these scorers use. Under the ACMG/AMP guidelines, computational predictions and evolutionary conservation count as supporting evidence for pathogenicity (criterion PP3) or benignity (BP4), and expert panels for several control genes specify in-silico tools in their gene-specific rules. Control-gene variants have been curated far more heavily than oocyte-gene variants, so conservation-based evidence has had more opportunity to shape their labels. This would inflate control-gene AUROC for any conservation-aware scorer, ESM-2 more than BLOSUM62. The review-status analysis shows the gap does not depend on well-reviewed variants, but 0-1 star submissions can also cite in-silico evidence, so it does not remove the problem. [After step 11: Deep mutational scanning labels for four control genes, which do not depend on curation, give ESM-2 AUROCs of [ ].]

### 5.3 New paragraph: negative-class annotation

> All variants removed by sequence validation were gnomAD negatives, which suggests that some gnomAD protein annotations refer to a non-canonical transcript. Validation checks only the wild-type residue, so a small number of retained negatives may be placed at the wrong position by chance. This adds noise to the negative class, most for genes with high drop rates (PTEN, TP53, PATL2, TUBB8).

### 5.4 Panel-size paragraph

**Before**
> Four oocyte genes clear the inclusion threshold and three of them carry fewer than ten pathogenic variants, so the oocyte panel mean rests on limited data even though the leave-one-out check indicates it is not driven by any single gene.

**After**
> Four oocyte genes clear the inclusion threshold and three of them carry fewer than ten pathogenic variants, so the oocyte panel mean rests on limited data even though the leave-one-out check indicates it is not driven by any single gene. With four oocyte and ten control genes, the permutation test can return no p-value below 1/1,001, and under the stricter frequency floor no value below 1/84.

---

## 6. Style

Words to remove or replace wherever they appear: "flagship" (use "most studied"), "clearest way to state", and any rhetorical framing such as "A model with 650 million parameters ... delivers approximately two hundredths". State the numbers and let them stand.
