# 06 — Motif discovery
## Reference and provenance

Parts of the motif discovery logic are inspired by or adapted from the
approach described in Sahu et al. (2022). The original reference
implementation is available at:

https://zenodo.org/records/5159644

## Purpose
Identify enriched sequence motifs associated with active inserts and/or TSS-associated elements.

## Inputs
- Signal sequences (active inserts, TSS-flanking sequences)
- Background sequences (input library or matched controls)
- Insert/activity tables (Step 04)

## Outputs
- Motif candidates (PWM/PPM)
- Motif rankings (support, enrichment)
- Comparison to known TF motifs (optional)

## Approach
1) Localmax-Motif: initial seed discovery  
2) Filter high-confidence seeds  
3) SeedExtender to extend seeds  
4) Rank/select extended motifs  
5) MEME suite de novo motif discovery  
6) PWM generation and reporting

## Code locations
- Localmax-Motif: github code
- SeedExtender:   github code
- MEME commands:  github code
- Spacek40:       github code
- Autoseed(compilation):       github code
## Notes / TODO
- Define final input sequence lengths (e.g., 150bp background, 150bp signal)
- Decide whether to use unique inserts only or weighted by counts/activity
- Sahu may heve used other methodes too
