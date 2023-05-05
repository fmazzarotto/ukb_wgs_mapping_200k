# ukb_wgs_mapping
This repo contains a map of the genomic coordinates for each WGS VCF block of the 200k WGS UK Biobank release, and the code used to create it.

The R script file "create_wgs_chunk_map.R" contains a single function which takes a table with three columns as input ("n_blocks_bases_per_chr.csv"). 
The input table is structured in a 1-row-per-chromosome basis, with column 1 containing the chromosome (with no chr prefix), column 2 containing 
the last UKB WGS block number for that chromosome (e.g. b4979 is the last block for chromosome 1, so the block number for 
chromosome 1 is 4979) and column 3 containing the last base number for the same chromosome (i.e. the chromosome length). 
Note that WGS blocks are numbered in a 0-based manner (i.e. chromosome 1 has 4980 blocks but the last block is b4979).

Each block contains 50kb of sequence, with the exception of regions 49100001-49150000 of chromosome 4 and 41850001-41900000 of 
chromosome 10, which have 5kb blocks. The function creates a table ("WGS_200k_block_map.tsv") which is a map of exact coordinates contained by 
each WGS block.
