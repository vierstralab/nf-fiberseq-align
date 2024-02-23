# Repo to align fiber-seq data

## main.nf
Align each flow cell separately + extract signal

## aggregation.nf
Merge individual alignments according to group column and divide them by chromosome


## installation notes. fibertools cannot be adequately installed with conda. Use cargo install --all-features --git https://github.com/fiberseq/fibertools-rs

