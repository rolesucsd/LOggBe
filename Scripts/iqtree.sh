#!/bin/zsh

source /home/roles/anaconda3/bin/activate

conda activate iqtree

iqtree -m GTR+I+G -s $1 -nt AUTO