# Knitwork
refactor of FragmentKnitwork

## Installation

```
git clone https://github.com/xchem/Knitwork
cd Knitwork
pip install --user -e .
```

## Configuration

The `configure` command can be used to set variables that the package will use when running commands:

e.g.

```
python -m knitwork configure GRAPH_LOCATION XXXXX
python -m knitwork configure GRAPH_USERNAME XXXXX
python -m knitwork configure GRAPH_PASSWORD XXXXX
```

## Running Fragment Knitwork

Run the following steps to generate merges from ligands in an SDF. 

## Fragmentation

To run the Fragment process which looks for subnodes and synthons for a given set of fragments/molecules in an SDF and groups them into pairs:

```
python -m knitwork fragment INPUT_SDF
```

This will generate pickled pandas dataframes, along with caches and other outputs in `fragment_output` by default:

- `molecules.pkl.gz`: pickled dataframe of input molecules
- `molecules.sdf`: SDF of input molecules
- `pairs.pkl.gz`: pickled dataframe of output pairs

For more options see:

```
python -m knitwork fragment --help
```

## Pure Knitting

To query the graph database for "pure" merges matching fragment pairs in the `fragment_output` folder by default:

```
python -m knitwork pure_merge
```

This will generate pickled pandas dataframes, along with caches and other outputs in `knitwork_output` by default:

- `pure_merges.pkl.gz`: pickled dataframe of merges
- `pure_merges.sdf`: SDF of merges

## Impure Knitting
