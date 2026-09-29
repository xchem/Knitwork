# Knitwork

[![latest](https://github.com/xchem/knitwork/actions/workflows/latest.yaml/badge.svg)](https://github.com/xchem/knitwork/actions/workflows/latest.yaml)
[![release](https://github.com/xchem/knitwork/actions/workflows/release.yaml/badge.svg)](https://github.com/xchem/knitwork/actions/workflows/release.yaml)

Refactor of FragmentKnitwork

## Installation

Knitwork is published to [PyPI](https://pypi.org/project/xchem-knitwork/) as `xchem-knitwork`
(the Python package is still imported as `knitwork`):

```
pip install xchem-knitwork
```

For development, install from a clone of the repository:

```
git clone https://github.com/xchem/Knitwork
cd Knitwork
pip install --user -e .
```

## Testing

Tests use [uv](https://docs.astral.sh/uv/) and the `dev` dependency group:

```
uv run --only-group dev pytest
```

## Releasing

Releases are made by pushing a semantic version tag (with no `v` prefix),
e.g. `1.2.0` or `1.2.0-rc.1`. The package version is taken from the tag.
The release workflow then publishes: -

- The Python package to PyPI (using PyPI Trusted Publishing)
- A Docker image to Docker Hub (as `xchem/knitwork:<tag>`)

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
python -m knitwork pure-merge
```

This will generate pickled pandas dataframes, along with caches and other outputs in `knitwork_output` by default:

- `pure_merges.pkl.gz`: pickled dataframe of merges
- `pure_merges.sdf`: SDF of merges

## Impure Knitting

To query the graph database for "impure" merges matching fragment pairs in the `fragment_output` folder by default:

```
python -m knitwork impure-merge
```

This will generate pickled pandas dataframes, along with caches and other outputs in `knitwork_output` by default:

- `impure_merges.pkl.gz`: pickled dataframe of merges
- `impure_merges.sdf`: SDF of merges
