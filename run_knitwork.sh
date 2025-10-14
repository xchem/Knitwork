
# fail anytime unset variables are used
set -u

# configure graph location
python -m knitwork.configure GRAPH_LOCATION $NEO4J_LOCATION
python -m knitwork.configure GRAPH_USERNAME $NEO4J_USERNAME
python -m knitwork.configure GRAPH_PASSWORD $NEO4J_PASSWORD

mkdir -p knitwork_input

# combine individual sdfs into a single input
python -m knitwork combine-inputs *.sdf --output knitwork_input/input.sdf

# run fragmentation on input.sdf
python -m knitwork fragment knitwork_input/input.sdf

# outputs: 
# - fragment_output/molecules.pkl.gz
# - fragment_output/molecules.sdf
# - fragment_output/pairs.pkl.gz

# run pure merging
python -m knitwork pure-merge

# outputs: 
# - knitwork_output/pure_merges.pkl.gz
# - knitwork_output/pure_merges.sdf

# run impure merging
python -m knitwork impure-merge

# outputs: 
# - knitwork_output/impure_merges.pkl.gz
# - knitwork_output/impure_merges.sdf