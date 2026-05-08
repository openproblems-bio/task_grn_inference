# SoundLife — GRN Inference & Evaluation Datasets

**Source:** `/vol/projects/jnourisa/hiara/datasets_old/sc/soundlife.h5ad`  
**Cohort:** SoundLife longitudinal flu-vaccination study  
**Cells:** 13,789,548 total PBMCs · 13,164 genes  
**Cell type used:** B cells (Major_CT == 'B') · ~1.2M B cells total  
**Donors per split:** 10 (shared across inference and evaluation)

---

## Temporal structure

| visitName | day | year | vaccinated |
|---|---|---|---|
| Flu Year 1 Day 0 | 0 | 1 | 0 (pre-vax) |
| Flu Year 1 Day 7 | 7 | 1 | 1 (post-vax) |
| Flu Year 1 Day 90 | 90 | 1 | 1 (post-vax) |
| Flu Year 2 Day 0 | 0 | 2 | 0 (pre-vax) |
| Flu Year 2 Day 7 | 7 | 2 | 1 (post-vax) |
| Flu Year 2 Day 90 | 90 | 2 | 1 (post-vax) |

Day 0 is confirmed pre-vaccination (vaccinated == 0 for all day 0 cells).

---

## Scenario 1a — `script_1.py`

**Cell type:** CD4T (subtypes: Tcm_Naive_CD4, Tem_Effector_CD4, Treg)  
**Question:** Does a GRN inferred from resting CD4T cells in year 1 generalise to the same donors one year later?

| Split | Condition | ~Cells |
|---|---|---|
| Inference | CD4T, day=0, year=1 | ~14K |
| Evaluation | CD4T, day=0, year=2 | ~13K |

**Outputs:**
- `resources/grn_benchmark/inference_data/soundlife_rna_s1.h5ad`
- `resources/grn_benchmark/evaluation_data/soundlife_sc_s1.h5ad`

---

## Scenario 2 — `script_2.py`

**Cell type:** B (subtypes: Naive_B, Memory_B, Plasma_B)  
**Question:** Does a GRN inferred from baseline + acute vaccination response (day 0 + day 7, year 1) predict the 90-day memory consolidation state across year 1 and year 2?

| Split | Condition | ~Cells |
|---|---|---|
| Inference | B, day=0+7, year=1 | ~28K |
| Evaluation | B, day=90, year=1+2 | ~23K |

**Outputs:**
- `resources/grn_benchmark/inference_data/soundlife_B_rna_s2.h5ad`
- `resources/grn_benchmark/evaluation_data/soundlife_B_sc_s2.h5ad`

---

## Output h5ad format

```
.X                  raw counts (float32, sparse CSR)
.layers['lognorm']  normalize_total(1e4) → log1p
.obs:
    donor_id        anonymised (donor_0 … donor_9)
    cell_type       'B'
    cell_type_orig  Sub_CT (Naive_B, Memory_B, Plasma_B, …)
    perturbation    pre_vaccination | post_vaccination_d7 | post_vaccination_d90
    day             '0' | '7' | '90'
    year            '1' | '2'
    visitName       original visit label
    donor_age       age at draw
    donor_sex       biological sex
.var.index          gene symbol (name = 'gene')
.uns:
    dataset_id      'soundlife'
    normalization_id 'lognorm'
    split           'inference' | 'evaluation'
    scenario        '1a' | '2'
```

---

## Run

```bash
# from task_grn_inference root
python src/process_data/main/soundlife/script_1.py
python src/process_data/main/soundlife/script_2.py
```
