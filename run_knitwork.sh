#!/bin/bash

set -u # fail anytime unset variables are used
set -e # exit on any failed commands

LIGANDS=$@
echo LIGANDS $LIGANDS

OUTPUT_DIR=$(dirname "$1")
echo OUTPUT_DIR $OUTPUT_DIR

CONFIG_PATH="knitwork_config.json"

mkdir -pv $OUTPUT_DIR/knitwork_input

### configure variables

python -m knitwork configure \
	GRAPH_LOCATION "bolt://"$NEO4J_LOCATION \
	--config-path=$CONFIG_PATH --silent

python -m knitwork configure \
	GRAPH_USERNAME $NEO4J_USERNAME \
	--config-path=$CONFIG_PATH --silent

python -m knitwork configure \
	GRAPH_PASSWORD $NEO4J_PASSWORD \
	--config-path=$CONFIG_PATH --silent

### combine individual sdfs into a single input

python -m knitwork combine-inputs \
	--root=data \
	$LIGANDS \
	$OUTPUT_DIR/knitwork_input/input.sdf

### run fragmentation on input.sdf

python -m knitwork fragment \
	$OUTPUT_DIR/knitwork_input/input.sdf \
	--config-path=$CONFIG_PATH \
	--output-dir=$OUTPUT_DIR/fragment_output

	# outputs: 
	# - fragment_output/molecules.pkl.gz
	# - fragment_output/molecules.sdf
	# - fragment_output/pairs.pkl.gz

### run pure merging

python -m knitwork pure-merge \
	--config-path=$CONFIG_PATH \
	--fragment-dir=$OUTPUT_DIR/fragment_output \
	--output-dir=$OUTPUT_DIR/knitwork_output

	# outputs: 
	# - knitwork_output/pure_merges.pkl.gz
	# - knitwork_output/pure_merges.sdf

### run impure merging

python -m knitwork impure-merge \
	--config-path=$CONFIG_PATH \
	--fragment-dir=$OUTPUT_DIR/fragment_output \
	--output-dir=$OUTPUT_DIR/knitwork_output

	# outputs: 
	# - knitwork_output/impure_merges.pkl.gz
	# - knitwork_output/impure_merges.sdf
