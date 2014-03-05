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


# read input
# drb1a1  drb1a2	dpb1a1	dpb1a2	dqb1a1	dqb1a2	dqa1a1	dqa1a2
# 04:02	07:01	72:01	132:01	02:02	03:02	02:01	03:00
# 03:01	07:01	01:01	02:01	02:01	03:03	02:01	05:01

setwd("~/Projects/em_haplotype")
#input_data = read.table("cleaned_calls.txt", sep="\t", header=T)
input_data = read.table("test.txt", sep="\t", header=T, stringsAsFactors=F)
marker_num = 4
sample_num = dim(input_data)[1]
head(input_data)

# Let G contain the unphased genotype data (input)
# fvec, be a vector of haplotype frequencies
# G is a  n x m x 2 array

master_haplotype = c()
G = array(data=NA, dim = c(sample_num, marker_num, 2))
for (n in 1:sample_num) {
  for (m in 1:marker_num) {
    G[n, m, 1] = input_data[n, 2*m-1]
    G[n, m, 2] = input_data[n, 2*m]
  }
  
  #ignore # N.a is a vector of length m, giving the number of unique alleles at each marker   
  
  # ENUMERATE the possible hapltoypes given the observed data
  #   to do so, go thru each individual and enumerate theirs,
  #   	adding each possible haplotype pair to the master haplotype
  #			list.  Then the dimenstion/length (LEN) of that list is the length of f
  
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.'))
  
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.'))
  master_haplotype = append(master_haplotype, paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.'))
}

# then we can initiate f to 1/ LEN
f = rep(1/length(master_haplotype), length(master_haplotype))
head(f)
length(f)

# length(unique(c( G[, 2 ,] ) ) )
# countunique <- function

## now we consider the E step

# H is a list, with an entry for each diploid patient/individual (length n)
H = list()

# each entry of H will in turn be a vector or list.  the number of such entries will be 2**(nh - 1), 
#where nh is the number of heterozygotes for that individual/entry

# for each person/patient
# 	foreach haplotype configuration
#		identify the two haplotypes (possibly via concatenation; ie. a string)
#		calculate the prob. of that configuration as the product of the two elements of f 
#   (corresponding to the two hap/strings defining the configuration) *TIMES* 2**(string1 != string2)
# bascially it means 2 choose 1
#      normalize by dividing by sum of probabilities

# permutations of configs
perm = matrix(c(1,1,1,1,1,1,1,2,1,1,2,1,1,1,2,2,1,2,1,1,1,2,1,2,1,2,2,1,1,2,2,2),nrow=8, ncol=4, byrow=T) 

# initialize entries of Exp_haps to zero.  Basically this will be just like your original "existence" dictionary;  master list of possible haplotypes given G
exp_haps = rep(0,dim(G)[1])
for (n in 1:dim(G)[1]) {
  
  H = append(H, n)
  
  # permute all possible haplotype configurations
  hconfig_sum = 0
  #H[[n]] = NULL
  for (p in 1:dim(perm)[1]) {
    haplotype1 = paste(c(G[n,1,perm[p,1]], G[n,2,perm[p,2]], G[n,3,perm[p,3]], G[n,4,perm[p,4]]), collapse='.')
    haplotype2 = paste(c(G[n,1,inv(perm[p,1])], G[n,2,inv(perm[p,2])], G[n,3,inv(perm[p,3])], G[n,4,inv(perm[p,4])]), collapse='.')
    f_haplotype1 = f[match(haplotype1, master_haplotype)]
    f_haplotype2 = f[match(haplotype2, master_haplotype)]
    
    if (haplotype1 == haplotype2) {
      f_haplotype_config = f_haplotype1 * f_haplotype2
    } else {
      f_haplotype_config = 2 * f_haplotype1 * f_haplotype2
    }
    hconfig_sum = hconfig_sum + f_haplotype_config
    
    #todo: set name 
    #af_haplotype_configs[paste(c(haplotype1,'/',haplotype2), collapse='')] = f_haplotype_config
    #H[[n]][p] = paste(c(haplotype1,'/',haplotype2,'-',f_haplotype_config), collapse='')
    H[[n]][p] = f_haplotype_config
    #names(H[[n]]) <- c(names(H[[n]]) , paste( sort( haplotype1, haplotype2), sep="/", collapse="" ) )  # same as line below
    #names(H[[n]])[ length(H[[n]]) ] <-  paste( sort( haplotype1, haplotype2), sep="/", collapse="" ) 
  }
  
  for (p in 1:dim(perm)[1]) {
    H[[n]][p]  <- H[[n]][p] / hconfig_sum
    
    # grab the two haplotypes associated with perm p (haplotype1, haplotype2)
    # add the expected haplotype count H[[n]][p] to Exp_haps for each of these 2 (!) haplotypes   
    # then after that, it is easy to do the M step 
  }
}

##  M step
#
# foreach element of f (ie, foreach haplotype possibly in our data set)
#    	count the number of EXPECTED times we see it  -- this is done in 2nd pound sign above -- simply grab the value associated with the key defined by this haplotype from Exp-haps
#			divide by (2*n)
#				that's our new estimate of f.

# stopping criteria

# when the likelihood changes less than delta #
# (we need to calculate this at each iteration)



# return compliment haplotype
inv <- function(x) {
  if (x==1) {
    return(2)
  } else {
    return(1)
  }
}


