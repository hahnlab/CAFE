seed 10

#specify data file and p-value threshold
load -i test3_families.txt -p 0.01

#the phylogenetic tree structure with branch lengths
tree (((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93)

#search for 2 parameter model
lambda -s -t (((2,2)1,(1,1)1)1,1)

#specify the global lambda to for generating simulated data
lambda -l 0.0017

rootdist -i fly.table

#generate 10 simulated data sets
genfamily rndtree/rnd -t 10

#estimate lambdas and compare likelihoods of global lambda and 2-parameter models
lhtest -d rndtree -l 0.0017 -t (((2,2)1,(1,1)1)1,1) -o lh2.out

