# HoGTIE (Horizontal Gene Transfer Identification Engine)
First thought to only occur widely in prokaryotic organisms, recent work has shown the prevalence of horizontal gene transfer (HGT) across the tree of life. HoGTIE uses a combination of k-mer-based and phylogenetic comparisons to identify genomic regions and sequences that contain HGT in a recipient species. 

### In development
HoGTIE is under active development. If you would like to check out or contribute to the code, hogtie can be installed locally via:

```
# coming soon...
# conda install hogtie -c conda-forge -c bioconda

# for now, install dependencies (toytree, numpy, pandas, scipy, sys, argparse) and
# do dev installation with pip

git clone https://github.com/cohen-r/hogtie
cd hogtie
pip install -e .

```

To test your installation, run:
```
hogtie --data sampledata/testdata.csv --tree sampledata/testtree.txt
```
The output should be:
```
Reading in data and tree...
Getting conditional likelihoods...
The conditional likelihoods at each node are:
[{0: 0.012983733049825068, 1: 0.013014488849433702}
 {0: 0.2175060189744712, 1: 0.24937595247663535}
 {0: 0.059509201574906405, 1: 0.05208798867475447}
 {0: 0.24999558915331158, 1: 0.24937595247663535}
 {0: 0.24999787050421876, 1: 0.24999815749065085}
 {0: 0.24533876067392912, 1: 0.2010535578911437}
 {0: 0.24999558915331163, 1: 0.24999588866552108}
 {0: 0.24533876067392907, 1: 0.24991613434302432}
 {0: 0.24991613434302426, 1: 0.24991613434302432} {0: 0, 1: 1}
 {0: 1, 1: 0} {0: 1, 1: 0} {0: 0, 1: 1} {0: 0, 1: 1} {0: 1, 1: 0}
 {0: 1, 1: 0} {0: 0, 1: 1} {0: 1, 1: 0} {0: 0, 1: 1}]
```