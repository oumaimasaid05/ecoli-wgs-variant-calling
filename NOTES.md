Day 1: setup
- Reference: E. coli K-12 MG1655, NC_000913.3
- Isolate: SRR28519659 (SRX24120536), ciprofloxacin-resistant colony 2, ZK126 + pRW320, PRJNA1095008
- Reads: Illumina MiniSeq, paired-end, 1,583,243 pairs, ~100x coverage, no subsampling
- Environment: environment.yml
- Day 3: aligned to MG1655 with bwa-mem2 (3,083,204 reads, ~62 s), sorted/indexed with samtools
- Day 4: variant calling with bcftools (ploidy 1): 749 raw, 44 filtered (38 SNPs, 6 indels); variant at 366125 flagged for strand bias
