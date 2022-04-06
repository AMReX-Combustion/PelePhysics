#!/usr/bin/env bash

# Prepare the QSS mechanism
mkdir qssaTools/mechanismFiles
cp skeletal.inp qssaTools/mechanismFiles
cp therm.dat qssaTools/mechanismFiles
cp tran.dat qssaTools/mechanismFiles
cp non_qssa_list.txt qssaTools/mechanismFiles
#cp non_qssa_list_save.txt qssaTools/mechanismFiles/non_qssa_list.txt

cd qssaTools

# Make graph
bash make_directed_graph.sh
echo "Made QSS Directed Graph"

cd ..

