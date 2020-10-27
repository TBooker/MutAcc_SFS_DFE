### Code to take the output of Remi's SimBit sims and to generate polyDFE configuration files
import argparse
import numpy as np
import pandas as pd
from collections import Counter

def createBoolArrayIsLocusNeutral(paramDict):
	nbLociPerCycle =  paramDict["nbCst"] + paramDict["nbGamma"] + paramDict["nbNtrl"] + paramDict["nbUnif"] 

	nbCycles = paramDict["nbLoci"] / nbLociPerCycle

#	isNeutral = np.array( [ ])

	isNeutral = np.repeat( np.concatenate( (np.repeat( paramDict["cstBool"], paramDict["nbCst"]),
				np.repeat( paramDict["gammaBool"], paramDict["nbGamma"]),
				np.repeat( paramDict["ntrlBool"], paramDict["nbNtrl"]),
				np.repeat( paramDict["unifBool"], paramDict["nbUnif"]) ), axis = None),
				nbCycles, axis = 0)

	if len(isNeutral) != paramDict["nbLoci"]:
		print("the number of loci does not match up with the expectation")
		return None
	else:
		return isNeutral
		
"""
R Code for the above function
  d = read.table(path, header=TRUE)
  
  nbLociperCycle = d$nbCst + d$nbGamma + d$nbNtrl + d$nbUnif
  nbCycles = d$nbLoci / nbLociperCycle
  
  isNeutral = rep(c(rep(d$cstBool, d$nbCst), rep(d$gammaBool, d$nbGamma), rep(d$ntrlBool, d$nbNtrl), rep(d$unifBool, d$nbUnif)), nbCycles)
  
  stopifnot(length(isNeutral) == d$nbLoci)
  return (isNeutral)
"""


def sampleSNPs(simbitOutput, neutralBool, sampleSize):
	for i in open(simbitOutput,"r"):
		if i.startswith("locus"):continue
		stringer = i.strip().split(" ")
		locus = int(stringer[0])
		freq = float(stringer[1])
		
		n_allele_copies = np.random.binomial(sampleSize, freq )
			
		if n_allele_copies == 0:
			continue
		else:
		#	print(n_allele_copies, freq, neutralBool[locus - 1])
			yield(n_allele_copies, neutralBool[locus])


def main():
	
## Start by defining the command line arguments, each one should be self-explanatory
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--simbit", 
		required = True,
		dest = "simbit",
		type = str, 
		help = "Give the simbit output file for this particular DFE")
	parser.add_argument("--param", 
		required = True,
		dest = "param",
		type = str, 
		help = "Give the parameter file for this particular DFE")
	parser.add_argument("--nInd", 
		required = True,
		dest = "nInd",
		type = int, 
		help = "What number of (haploid) individuals do you want to summarise?")
	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type = str, 
		help = "Give the name of the output file")
		
	args = parser.parse_args()
	
	param_DF =  pd.read_csv(args.param, sep = " ") 

	param_dict = param_DF.to_dict( "records")[0]  	

	isNeutral = createBoolArrayIsLocusNeutral( param_dict )

	synonymous_sites = isNeutral.sum()
	nonsynonymous_sites = isNeutral.shape[0] - synonymous_sites

## Now iterate through the SNP output file spit out the SNP freqs

	synonymous_counter = Counter()
	nonsynonymous_counter = Counter()

	for i in range(args.nInd + 1):
		synonymous_counter[i] = 0
		nonsynonymous_counter[i] = 0

	for s in sampleSNPs(args.simbit, isNeutral, args.nInd):
		if s[1] == 0.0: ## If s[1] == 0 this SNP is nonsynonymous
			nonsynonymous_counter[s[0]] +=1
		elif s[1] == 1.0:  ## If s[1] == 0 this SNP is synonymous
			synonymous_counter[s[0]] +=1
		
	syn_nonancestral = sum( synonymous_counter.values()) 
	nonsyn_nonancestral = sum( nonsynonymous_counter.values() )
	
	synonymous_counter[0] = synonymous_sites - syn_nonancestral
	nonsynonymous_counter[0] = nonsynonymous_sites - nonsyn_nonancestral
	
	SFS_df = pd.DataFrame.from_dict([ synonymous_counter, nonsynonymous_counter] ).transpose() 
	SFS_df.reset_index(inplace = True)
	SFS_df = SFS_df.rename( columns = {"index":"alleles", 0:"synonymous", 1:"nonsynonymous"})

	SFS_df.sort_values(["alleles"]).to_csv( args.output, index = False)


main()
