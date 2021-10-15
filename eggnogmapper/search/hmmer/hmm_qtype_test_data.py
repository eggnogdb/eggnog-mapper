##
# CPCantalapiedra 2021

test_hmm = """
HMMER3/f [3.1b2 | February 2015]
NAME  COG0012.faa.final_tree.fa
LENG  10
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
DATE  Fri Aug  3 19:28:36 2018
NSEQ  4522
EFFN  4.170202
CKSUM 2788465937
STATS LOCAL MSV      -12.2126  0.69661
STATS LOCAL VITERBI  -13.1336  0.69661
STATS LOCAL FORWARD   -6.4930  0.69661
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.51306  4.33786  2.93844  2.59829  3.23556  2.91764  3.69806  2.83251  2.59931  2.46112  3.66933  3.04067  3.48751  3.06333  2.89831  2.71431  2.83756  2.65132  4.65738  3.46485
          2.68513  4.41795  2.77550  2.73197  3.46434  2.40560  3.72429  3.29316  2.67903  2.69357  4.24254  2.90434  2.73762  3.18206  2.89505  2.37899  2.77572  2.98544  4.58540  3.61427
          1.15084  2.23638  0.55029  3.71384  0.02469  0.00000        *
      1   2.68937  4.25841  3.94611  3.41262  3.39439  3.71285  4.09869  2.64333  3.25914  2.42420  1.47597  3.66097  4.06968  3.56631  3.45847  2.95721  2.85866  2.11796  4.90994  3.69104    126 m - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01791  4.42731  5.14966  0.61958  0.77255  1.16884  0.37212
      2   1.99709  4.73712  3.21241  2.65685  3.97328  2.44917  3.75624  3.35257  2.61920  2.99654  3.11310  3.04381  3.74309  2.96868  3.07975  2.11439  2.90197  3.06289  5.30787  4.00367    127 a - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01677  4.49266  5.21500  0.61958  0.77255  1.00349  0.45665
      3   2.98280  4.35623  4.80266  4.13824  2.47496  4.19294  4.51141  2.38493  4.00936  1.16926  1.96464  4.30877  4.51761  4.18939  4.10745  3.50375  3.21203  2.47025  4.86770  3.79496    128 l - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01349  4.70806  5.43041  0.61958  0.77255  1.11135  0.39916
      4   2.71960  5.25207  2.98103  2.37043  4.60202  3.45558  3.69218  4.07289  1.64831  3.56423  4.26072  2.66771  3.94258  2.43515  2.73206  2.42250  2.38468  3.65083  5.69849  4.21908    129 k - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01347  4.70987  5.43221  0.61958  0.77255  1.11135  0.39916
      5   2.10517  2.24730  4.70335  4.10126  3.29896  4.02021  4.35980  1.84453  3.91602  2.09737  2.75921  4.05672  4.38038  4.07207  3.96895  3.31760  2.99406  2.03557  4.88488  3.69568    130 i - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01320  4.73000  5.45235  0.61958  0.77255  1.11135  0.39916
      6   2.28072  4.29893  4.34000  4.09789  4.97337  0.53875  4.99787  4.38122  4.04976  3.95508  4.91350  3.98340  4.18203  4.32598  4.02830  2.79913  3.15752  3.67911  6.32659  5.17976    131 G - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01320  4.73000  5.45235  0.61958  0.77255  1.11135  0.39916
      7   3.25840  4.61370  4.92706  3.77016  3.46589  4.50211  4.85936  0.80595  4.20572  1.89747  3.45746  4.53357  4.80672  4.28707  4.34038  3.83250  3.49477  2.19928  5.37473  4.19960    132 i - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01472  4.73025  5.15268  0.61958  0.77255  1.11135  0.39916
      8   3.22508  4.55192  5.32528  4.79429  3.74838  4.79912  5.30799  1.65421  4.68133  2.39109  3.73677  4.91868  5.08319  4.88926  4.80400  4.16543  3.52153  0.70717  5.76062  4.30935    133 v - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01319  4.73081  5.45316  0.61958  0.77255  1.11307  0.39831
      9   2.93961  4.89785  3.82919  3.77726  4.95651  0.39957  4.63387  4.68860  3.98064  4.32347  4.94405  3.95042  4.37067  4.29848  4.22942  3.13906  3.49582  4.05851  6.29839  5.01721    134 G - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01316  4.73273  5.45508  0.61958  0.77255  1.06255  0.42399
     10   3.18374  4.48400  4.89614  4.18336  2.80495  4.36955  4.65328  2.64311  4.05799  0.70915  2.99013  4.43698  4.61782  4.27264  4.22131  3.62110  3.42556  2.60780  5.10292  3.83008    135 l - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01598  4.74807  4.93623  0.61958  0.77255  1.09067  0.40946
//
"""

## END