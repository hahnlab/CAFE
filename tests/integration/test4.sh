seed 10

tree ((((cat:68,horse:68):4,cow:73):20,(((((chimp:4,human:4):6,orang:11):2,gibbon:13):7,(macaque:4,baboon:4):16):16,marmoset:36):57):38,(rat:36,mouse:36):96)

load -i test4_families.txt

errormodel -all -model errormodel.txt

lambda -l 0.01 0.005 -t ((((1,1)1,1)1,(((((2,2)2,2)2,2)2,(1,1)1)1,1)1)1,(1,1)1) -score

