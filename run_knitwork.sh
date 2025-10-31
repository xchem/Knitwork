#!/bin/bash

set -u # fail anytime unset variables are used
set -e # exit on any failed commands

echo CMD $0 $@
echo PWD $(pwd)
echo LS $(ls -la)

CONFIG_PATH="knitwork_input/knitwork_config.json"

mkdir -pv knitwork_input

# configure graph location
python -m knitwork configure GRAPH_LOCATION $NEO4J_LOCATION --config-path=$CONFIG_PATH
python -m knitwork configure GRAPH_USERNAME $NEO4J_USERNAME --config-path=$CONFIG_PATH
python -m knitwork configure GRAPH_PASSWORD $NEO4J_PASSWORD --config-path=$CONFIG_PATH

# combine individual sdfs into a single input
python -m knitwork combine-inputs data/*.sdf knitwork_input/input.sdf

# run fragmentation on input.sdf
python -m knitwork fragment knitwork_input/input.sdf --config-path=$CONFIG_PATH

# outputs: 
# - fragment_output/molecules.pkl.gz
# - fragment_output/molecules.sdf
# - fragment_output/pairs.pkl.gz

# run pure merging
python -m knitwork pure-merge --config-path=$CONFIG_PATH

# outputs: 
# - knitwork_output/pure_merges.pkl.gz
# - knitwork_output/pure_merges.sdf

# run impure merging
python -m knitwork impure-merge --config-path=$CONFIG_PATH

# outputs: 
# - knitwork_output/impure_merges.pkl.gz
# - knitwork_output/impure_merges.sdf