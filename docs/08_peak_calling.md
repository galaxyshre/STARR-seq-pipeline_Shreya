# 08 — Genome-wide activity peak calling (optional)

## Purpose
Identify contiguous genomic regions with elevated regulatory activity by aggregating STARR-seq RNA/DNA enrichment into genome-wide activity tracks and defining active regions ("peaks"), following the conceptual framework described in **:contentReference[oaicite:0]{index=0}**.

This step is **optional** and applies **only to genome-aligned STARR-seq libraries**
(e.g. genomic promoters or enhancers).  
Random sequence libraries are excluded.

---

## Scope and rationale

STARR-seq measures transcriptional output rather than protein binding.  
Therefore, peak calling is not intrinsic to activity quantification and is **not required** for:
- insert-level activity estimation
- motif discovery
- promoter grammar modeling
- sequence-to-activity learning

Peak calling is included here as an **auxiliary analysis** for:
- genome browser visualization
- comparison with external genomic datasets
- annotation of active regulatory regions

Insert-level activity (Step 04) remains the primary quantitative signal.

---

## Inputs
- Genome-aligned STARR-seq data:
  - RNA BAM files (signal)
  - DNA input BAM files (background)
- Genome reference (e.g. hg38)
- Replicate metadata
- (Optional) blacklist regions or mappability filters

---

## Outputs
- Genome-wide activity tracks:
  - normalized RNA coverage
  - normalized DNA coverage
  - activity track (log2 RNA/DNA)
- Peak files (BED format), containing:
  - peak_id
  - genomic coordinates
  - peak length
  - summary activity score (mean or max)
- Replicate-level and merged peak sets

---

## Method (conceptual framework)

Peak calling is performed **after insert-level activity quantification** and consists of the following steps:

1. **Generate coverage tracks**  
   Compute genome-wide RNA and DNA coverage from aligned reads and normalize by library size.

2. **Compute activity signal**  
   For each base or genomic bin:
   \[
   \text{activity} = \log_2 \left( \frac{\text{RNA} + \epsilon}{\text{DNA} + \epsilon} \right)
   \]
   where \(\epsilon\) is a small pseudocount to stabilize low-coverage regions.

3. **Signal smoothing**  
   Smooth the activity track using sliding windows or kernel-based smoothing to reduce local noise and emphasize contiguous regulatory domains.

4. **Thresholding**  
   Define active regions using empirical thresholds, such as:
   - percentile-based cutoffs (e.g. top X% of activity)
   - Z-score relative to genomic background
   - comparison to shuffled or inactive regions

5. **Peak definition**  
   Merge adjacent above-threshold bins into contiguous regions to define peaks.

6. **Replicate handling**  
   Assess reproducibility across replicates and generate merged peak sets using overlap or consensus-based approaches.

Exact parameters (bin size, smoothing window, thresholds) may be tuned per dataset.

---

## QC checks / expected results
- Activity tracks show structured enrichment rather than sparse spikes
- Peaks cluster near known regulatory regions (e.g. promoters or enhancers)
- Peak widths are consistent with promoter/enhancer architecture
- Reproducible peak calls across biological replicates
- Enrichment near known TSSs and active chromatin marks (where applicable)

---

## Relationship to other pipeline steps
- Depends on: **Step 04 — Activity quantification**
- Independent of:
  - motif discovery
  - TSS inference
  - sequence-based modeling
- Peak calling does **not** feed back into insert-level activity estimates

---

## Notes
- Standard ChIP-seq peak callers (e.g. MACS2) are not directly applicable to STARR-seq data.
- STARR-seq peak calling reflects **transcriptional activity**, not factor binding.
- Random sequence libraries cannot be subjected to this analysis due to lack of genomic coordinates.

---

## Reference
This analysis follows the conceptual approach described in:  
**:contentReference[oaicite:1]{index=1}**

