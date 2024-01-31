#!/bin/bash

i=$((RANDOM % ${1}))
echo "sleeping ${i}"
sleep ${i} 

echo "installing"
pip install --user h5py
pip install --user numpy
pip install --user scipy
