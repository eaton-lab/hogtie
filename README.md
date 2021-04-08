# HoGTIE (**Ho**rizontal **G**ene **T**ransfer **I**dentification **E**ngine)
First thought to only occur widely in prokaryotic organisms, recent work has shown the prevalence of horizontal gene transfer (HGT) across the tree of life. HoGTIE uses a Discrete Markov Model to reconstruct ancestral character states of sequence motifs (*i.e.* kmers, SNPs, or transcripts) across a species tree, identifying sequence motifs that are likely a result of introgression and mapping them along a genome-of-interest. HoGTIE employs a flexible, comparative GWAS method that can be used in the absence of well-resolved genome assemblies.

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
hogtie --data sampledata/testdata.csv --tree sampledata/testmatrix.txt
```