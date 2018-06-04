#!~/bin/cafe
seed 10
load -i test2_families.txt -p 0.05 -max_size 20
tree (((Chimp:6,Human:6):81,(Mouse:17,Rat:17):70):6,Dog:93)
lambda -s
report test2

