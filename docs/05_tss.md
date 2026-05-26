# 05 — TSS identification about RPRE library
## Reference implementation (Sahu et al.)

TSSs were identified by mapping the TSO library back to the input library, using a strategy similar to that described by Sahu et al.

## Purpose
Identify transcription start sites (TSSs) captured by template-switch (TS) libraries and associate TSS positions with inserts.

## Inputs
- TSL library FASTQs (R1)
- Random library 
- The code written by Dan, Shreya, Farzin and Fei. 

## Outputs
- TSS calls 
- Summary stats: number of TSS reads, adapter-containing fraction, unique TSS counts

## Notes / TODO
- Specify read structure assumptions (TS adapter, primers, UMI if any)
- Define final output formats expected by downstream modelling

