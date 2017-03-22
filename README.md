# CAFE

Software for **C**omputational **A**nalysis of gene **F**amily **E**volution

The purpose of CAFE is to analyze changes in gene family size in a way that 
accounts for phylogenetic history and provides a statistical foundation for 
evolutionary inferences. The program uses a birth and death process to model gene 
gain and loss across a user-specified phylogenetic tree. The distribution of family 
sizes generated under this model can provide a basis for assessing the significance 
of the observed family size differences among taxa.

CAFE v4.0 is the first in a regular series of releases to the CAFE application. The 
manual and various tutorials may be viewed on the website. This document describes how to 
download and use CAFE v4.0. 

# Use

The necessary inputs for CAFE v4.0 are:
1.  a data file containing gene family sizes for the taxa included in the 
phylogenetic tree
2.  a Newick formatted phylogenetic tree, including branch lengths

From the inputs above, CAFE v4.0 will compute:
1.  the maximum likelihood value of the birth & death parameter, λ (or of 
separate birth and death parameters (λ and μ, respectively), over the whole 
tree or for user-specified subsets of branches in the tree
1.  ancestral states of gene family sizes for each node in the phylogenetic tree
1.  p-values for each gene family describing the likelihood of the observed sizes 
given average rates of gain and loss
1.  average gene family expansion along each branch in the tree
1.  numbers of gene families with expansions, contractions, or no change
along each branch in the tree

# Install

Run "make" from the home directory. If the readline library is available on the system, running
"make USE_READLINE=1" will allow CAFE to support a command line history.

# History

CAFE v3.0 was a major update to CAFE v2.1. Major updates in 3.0 included: 1) the ability to correct 
for genome assembly and annotation error when analyzing gene family evolution 
using the **errormodel** command. 2) The ability to estimate separate birth (λ) and 
death (μ) rates using the **lambdamu** command. 3) The ability to estimate error in 
an input data set with iterative use of the errormodel command using the 
accompanying python script **caferror.py**. This version also included the addition of the 
**rootdist** command to give the user more control over simulations.

