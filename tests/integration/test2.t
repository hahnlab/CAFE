  $ cafe -v
  Version: 4.1, built at Apr 11 2018

  $ cd $TESTDIR && cafe test2.sh
  -----------------------------------------------------------
  Family information: test2_families.txt
  Log: stdout
  The number of families is 4
  Root Family size : 1 ~ 30
  Family size : 0 ~ 56
  P-value: 0.05
  Num of Threads: 1
  Num of Random: 1000
  (((Chimp:6,Human:6):81,(Mouse:17,Rat:17):70):6,Dog:93)
  Empirical Prior Estimation Result: (33 iterations)
  Poisson lambda: 1.700000 & Score: 34.151494
  Lambda : 0.00656913889832 & Score: -31.227238
  .Lambda : 0.00689759584324 & Score: -31.629977
  .Lambda : 0.00624068195341 & Score: -30.814590
  .Lambda : 0.00591222500849 & Score: -30.391758
  .Lambda : 0.00525531111866 & Score: -29.514961
  .Lambda : 0.00459839722883 & Score: -28.597658
  .Lambda : 0.00328456944916 & Score: -26.679683
  .Lambda : 0.00197074166950 & Score: -24.946755
  .Lambda : -0.00065691388983 & Score: -inf
  .Lambda : 0.00065691388983 & Score: -25.619505
  .Lambda : 0.00328456944916 & Score: -26.679683
  .Lambda : 0.00262765555933 & Score: -25.743380
  .Lambda : 0.00131382777966 & Score: -24.575773
  .Lambda : 0.00065691388983 & Score: -25.619505
  .Lambda : 0.00098537083475 & Score: -24.784846
  .Lambda : 0.00164228472458 & Score: -24.675955
  .Lambda : 0.00098537083475 & Score: -24.784846
  .Lambda : 0.00114959930721 & Score: -24.628409
  .Lambda : 0.00147805625212 & Score: -24.598124
  .Lambda : 0.00114959930721 & Score: -24.628409
  .Lambda : 0.00123171354344 & Score: -24.591192
  .Lambda : 0.00139594201589 & Score: -24.578960
  .Lambda : 0.00123171354344 & Score: -24.591192
  .Lambda : 0.00127277066155 & Score: -24.580972
  .Lambda : 0.00135488489778 & Score: -24.575217
  .Lambda : 0.00139594201589 & Score: -24.578960
  .Lambda : 0.00137541345684 & Score: -24.576571
  .Lambda : 0.00133435633872 & Score: -24.574936
  .Lambda : 0.00131382777966 & Score: -24.575773
  .Lambda : 0.00132409205919 & Score: -24.575212
  .Lambda : 0.00134462061825 & Score: -24.574940
  .Lambda : 0.00132409205919 & Score: -24.575212
  .Lambda : 0.00132922419896 & Score: -24.575039
  .Lambda : 0.00133948847849 & Score: -24.574903
  .Lambda : 0.00134462061825 & Score: -24.574940
  .Lambda : 0.00134205454837 & Score: -24.574913
  .Lambda : 0.00133692240860 & Score: -24.574911
  .Lambda : 0.00134205454837 & Score: -24.574913
  .Lambda : 0.00134077151343 & Score: -24.574906
  .Lambda : 0.00133820544355 & Score: -24.574905
  .Lambda : 0.00134077151343 & Score: -24.574906
  .Lambda : 0.00134012999596 & Score: -24.574904
  .
  Lambda Search Result: 21
  Lambda : 0.00133948847849 & Score: 24.574903
  DONE: Lambda Search or setting, for command:
  lambda -s 
  Running Viterbi algorithm....
  Report Done

  $ cat test2.cafe
  Tree:(((Chimp:6,Human:6):81,(Mouse:17,Rat:17):70):6,Dog:93)
  (((Chimp:6,Human:6):81,(Mouse:17,Rat:17):70):6,Dog:93)Lambda:\t0.00133949 (esc)
  # IDs of nodes:(((Chimp<0>,Human<2>)<1>,(Mouse<4>,Rat<6>)<5>)<3>,Dog<8>)<7>
  # Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): (0,2) (1,5) (4,6) (3,8) 
  # Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (0, 1, 2, 3, 4, 5, 6, 7, 8)
  Average Expansion:\t(0.25,-0.25)\t(0.25,0.25)\t(0,0.25)\t(0,-0.25) (esc)
  Expansion :\t(1,0)\t(1,1)\t(0,1)\t(0,0) (esc)
  nRemain :\t(3,3)\t(3,3)\t(4,3)\t(4,3) (esc)
  nDecrease :\t(0,1)\t(0,0)\t(0,0)\t(0,1) (esc)
  'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'\t'cut P-value'\t'Likelihood Ratio' (esc)
  ENSF00000000007\t(((Chimp_4:6,Human_4:6)_4:81,(Mouse_3:17,Rat_3:17)_3:70)_3:6,Dog_3:93)_3\t0.6955\t((-,-),(-,-),(-,-),(-,-))\t (esc)
  ENSF00000000014\t(((Chimp_5:6,Human_3:6)_4:81,(Mouse_5:17,Rat_6:17)_5:70)_4:6,Dog_3:93)_4\t0.003\t((0.0159398,0.0461031),(0.749143,0.165864),(0.596659,0.0553873),(0.530652,0.419912))\t (esc)
  ENSF00000000015\t(((Chimp_1:6,Human_1:6)_1:81,(Mouse_1:17,Rat_1:17)_1:70)_1:6,Dog_1:93)_1\t0.7745\t((-,-),(-,-),(-,-),(-,-))\t (esc)
  ENSF00000000029\t(((Chimp_2:6,Human_2:6)_2:81,(Mouse_2:17,Rat_2:17)_2:70)_2:6,Dog_2:93)_2\t0.8985\t((-,-),(-,-),(-,-),(-,-))\t (esc)

  $ rm test2.cafe
 
