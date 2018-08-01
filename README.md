## Goal
* Previously have conducted single marker associations with scleroderma
* Want to test hypothesis that haplotypes may confer risk not attributable to additive effects of individual alleles
* Focuse on determining haplotypes of 4 HLA class II alleles (drb1, dpb1, dqb1, d1a1) using EM algorithm

## EM initialization
* Given 4 loci, each loci has 2 alleles (locus1 has A,a, locus2 has B,b, locus3 has C,c, locus4 has D,d)
* 739 individuals, N
* Up to 16 possible haplotypes per individual, h = {ABCD, ABCd, ABcD, …, abcd}
* 8 possible haplotype pairs per inidividual, { {ABCD,abcd}, {ABCd,abcD} … }
* Initial haplotype probabilities, f = (PABCD, PABCd, PABcD, PABcd, … ,Pabcd)
  + Random f
  + Equally likely f

## EM process
* Each individual i, has at most 8 haplotype pairs, each pair consists of haplotype j, and haplotype k. Conditional probability of a haplotype pair is

P_{i}(haplotype_pair = j,k | genotype_{AaBbCcDd}) proportion to F_{h_j} * f_{h_k} * 2^(h_j!=h_k)

* Then, the estimated frequency of a haplotype hj is the sum of its expected occurrences from all the individuals


