# EM for multilocus / multi-allelic haplotypes (eg HLA)
#
# Baisc steps:

# 1) Initiate with some haplotype allele frequencies, f

				# E step #
# 2) Given f, estimate the possible haplotype configurations
#  		for each individual, given their observed genotypes,
#			calculate the probability of each haplotype configuration

				# M step #
# 3) Estimate f', by considering the weighted average, based on the 
#			probabilities of the configurations.  
#  		Essentially, count the expected number of haplotypes (of each type) and use that to estimate f'

# 4)  set f = f', and go to (2)

#################################

# Let G contain the unphased genotype data (input)
# fvec, be a vector of haplotype frequencies
setwd("~/Projects/em_haplotype")
a = read.table("cleaned_calls.txt", sep="\t", header=T)
head(a)

all_hapotypes = 
# G is a  n x m x 2 array

#ignore # N.a is a vector of length m, giving the number of unique alleles at each marker   

# ENUMERATE the possible hapltoypes given the observed data
#   to do so, go thru each individual and enumerate theirs,
# 		adding each possible haplotype pair to the master haplotype
#			list.  Then the dimenstion/length (LEN) of that list is the length of f

# then we can initiate f to 1/ LEN

length(unique(c( G[, 2 ,] ) ) )
# countunique <- function

## now we consider the E step

# H is a list, with an entry for each diploid patient/individual (length n)

# each entry of H will in turn be a vector or list.  the number of such entries will be 2**(nh - 1), where nh is the number of heterozygotes for that individual/entry

# for each person/patient
# 	foreach haplotype configuration
#		identify the two haplotypes (possibly via concatenation; ie. a string)
#		calculate the prob. of that configuration as the product of the two elements of f (corresponding to the two hap/strings defining the configuration) *TIMES* 2**(string1 != string2)
#      normalize by dividing by sum of probabilities


##  M step

# foreach element of f (ie, foreach haplotype possibly in our data set)
#  		count the number of EXPECTED times we see it 
#			divide by (2*n)
#				that's our new estimate of f.

# set f' to f  #  f <- f'  and iterate.

