#!/bin/zsh

source /Users/reneeoles/opt/anaconda3/etc/profile.d/conda.sh

conda activate iqtree

iqtree -m GTR+I+G -s $1 -nt AUTO