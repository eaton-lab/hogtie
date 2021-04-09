# HoGTIE (**Ho**rizontal **G**ene **T**ransfer **I**dentification **E**ngine)
First thought to only occur widely in prokaryotic organisms, recent work has shown the prevalence of horizontal gene transfer (HGT) across the tree of life. HoGTIE uses a Discrete Markov Model to reconstruct ancestral character states of sequence motifs (*i.e.* kmers, SNPs, or transcripts) across a species tree, identifying sequence motifs that are likely a result of introgression and identifying regions that contain introgression along a linear genome. HoGTIE employs a flexible, comparative phylogenomic method that can be used in the absence of well-resolved genome assemblies.

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

If interested in using HoGTIE on a CLI, test your installation by running:
```
hogtie --tree sampledata/testtree.txt --data sampledata/testmatrix.csv --model ARD
```
HoGTIE will read in the testdata and write the results to a file in an output directory.

To use HoGTIE on the API (currently includes visualization of likelihood scores and ancestral character states along a tree), open up the API [working example in the notebooks folder](notebooks/working_example).