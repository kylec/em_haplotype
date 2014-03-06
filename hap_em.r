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

# return compliment haplotype
inv <- function(x) {
  if (x==1) {
    return(2)
  } else {
    return(1)
  }
}

# read input
# drb1a1  drb1a2	dpb1a1	dpb1a2	dqb1a1	dqb1a2	dqa1a1	dqa1a2
# 04:02	07:01	72:01	132:01	02:02	03:02	02:01	03:00
# 03:01	07:01	01:01	02:01	02:01	03:03	02:01	05:01

setwd("~/Projects/em_haplotype")
input_data = read.table("cleaned_calls.txt", sep="\t", header=T, stringsAsFactors=F, colClasses=c(rep("character",9)))
#input_data = read.table("test.txt", sep="\t", header=T, stringsAsFactors=F, colClasses=c(rep("character",9)))
marker_num = 4
sample_num = dim(input_data)[1]
rownames(input_data) = input_data[,1]
input_data = input_data[,c(2:9)]
head(input_data)

# Let G contain the unphased genotype data (input)
# f be a vector of haplotype frequencies
# G is a  n x m x 2 array
# exp_haps is a vector of expected number for a haplotype
#
# ENUMERATE the possible hapltoypes given the observed data
#   to do so, go thru each individual and enumerate theirs,
#     adding each possible haplotype pair to the master haplotype
#			list.  Then the dimenstion/length (LEN) of that list is the length of f

exp_haps = c()
G = array(data=NA, dim = c(sample_num, marker_num, 2))
f = c()

for (n in 1:sample_num) {
  for (m in 1:marker_num) {
    G[n, m, 1] = input_data[n, 2*m-1]
    G[n, m, 2] = input_data[n, 2*m]
  }
  
  f[paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.')] = 1
  f[paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.')] = 1
  
  exp_haps[paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,1], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,1], G[n,3,1], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,1], G[n,3,2], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,2], G[n,3,1], G[n,4,2]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,1]), collapse='.')] = 0
  exp_haps[paste(c(G[n,1,2], G[n,2,2], G[n,3,2], G[n,4,2]), collapse='.')] = 0
}

# initiate f to 1/ LEN
#x <- (runif(length(f)) + 1 ) 
#f = f*x/sum(x)
f = f/length(f)

## E step
#
# H is a list, with an entry for each diploid patient/individual (length n)
# each entry of H will in turn be a vector or list.  the number of such entries will be 2**(nh - 1), 
# where nh is the number of heterozygotes for that individual/entry

# for each person/patient
# 	foreach haplotype configuration
#		identify the two haplotypes (possibly via concatenation; ie. a string)
#		calculate the prob. of that configuration as f = p1 * p2 * 2 ** (hap1 != hap2)
#   normalize by dividing by sum of probabilities

# permutations of configs
perm = matrix(c(1,1,1,1,1,1,1,2,1,1,2,1,1,1,2,2,1,2,1,1,1,2,1,2,1,2,2,1,1,2,2,2),nrow=8, ncol=4, byrow=T) 

logdlikevec <- NULL
num_iter = 10
total_iter = 0

# do multiple em_step
for ( j in 1:1) {
# loop for em step
for (i in 1:num_iter) {
    # H is a list of samples, each sample has a list of haplotype and probability
    H = list()
    dlike = 0
    for (n in 1:dim(G)[1]) {
  
      H = append(H, n)
  
      # permute all possible haplotype configurations
      hconfig_sum = 0
      #H[[n]] = NULL
      for (p in 1:dim(perm)[1]) {
        haplotype1 = paste(c(G[n,1,perm[p,1]], G[n,2,perm[p,2]], G[n,3,perm[p,3]], G[n,4,perm[p,4]]), collapse='.')
        haplotype2 = paste(c(G[n,1,inv(perm[p,1])], G[n,2,inv(perm[p,2])], G[n,3,inv(perm[p,3])], G[n,4,inv(perm[p,4])]), collapse='.')
        f_haplotype1 = f[haplotype1]
        f_haplotype2 = f[haplotype2]
    
        # calculate haplotype config probability
        if (haplotype1 == haplotype2) {
          f_haplotype_config = f_haplotype1 * f_haplotype2
        } else {
          f_haplotype_config = 2 * f_haplotype1 * f_haplotype2
        }
        hconfig_sum = hconfig_sum + f_haplotype_config
    
        # haplotype config string is "key", probability is value
        H[[n]][p] = f_haplotype_config
        names(H[[n]])[p] = paste( sort( c(haplotype1, haplotype2)), collapse="/" )
    
      }
  
      dlike_n = 0
      # normalize config probability for a sample by sum
      for (p in 1:dim(perm)[1]) {
        dlike_n = dlike_n + H[[n]][p]
        H[[n]][p]  <- H[[n]][p] / hconfig_sum
    
        # grab the two haplotypes associated with perm p (haplotype1, haplotype2)
        # add the expected haplotype count H[[n]][p] to Exp_haps for each of these 2 (!) haplotypes   
        # then after that, it is easy to do the M step
        haplotype1 = unlist(strsplit(names(H[[n]])[p], "/"))[1]
        haplotype2 = unlist(strsplit(names(H[[n]])[p], "/"))[2]

        exp_haps[haplotype1] = exp_haps[haplotype1] + H[[n]][p]
        exp_haps[haplotype2] = exp_haps[haplotype2] + H[[n]][p]
      }
  
      dlike = dlike + log(dlike_n)
    }
    
    logdlikevec <- c(logdlikevec, dlike)
    ## M step
    #
    # foreach element of f (ie, foreach haplotype possibly in our data set)
    #    	count the number of EXPECTED times we see it  -- this is done in 2nd pound sign above
    #     simply grab the value associated with the key defined by this haplotype from Exp-haps
    #			divide by (2*n), new estimate of f

    for (i in 1:length(exp_haps)) {
      f_new = exp_haps[i]/(2*sample_num)
      f[names(exp_haps[i])] = f_new
    }
    head(f)
    head(exp_haps)
}

total_iter = total_iter + num_iter

#png(paste(total_iter, "iter.png", sep="_"))
plot(logdlikevec)
#dev.off()

# stopping criteria
# when the likelihood changes less than delta #
# (we need to calculate this at each iteration)

# select most likely haplotype config
output=c()
for (i in 1:length(H)) {
  max_hap_config = unlist(strsplit(names(which.max(H[[i]])),"/"))
  output = append(output, paste(rownames(input_data[i,]), max_hap_config[1], max_hap_config[2], sep="\t"))
}
write(output,file=paste(total_iter, "iter_haplotype_configs.txt", sep="_"))

}


