#!/bin/bash

echo $1 $2 $3 $4

biom convert -i $1 -o $2 --table-type='OTU table' --to-hdf5

faithpd -i $2 -o $3 -t $4