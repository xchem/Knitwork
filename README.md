# Knitwork
refactor of FragmentKnitwork

## Configuration

The `configure` command can be used to set variables that the package will use when running commands:

e.g.

```
python -m knitwork configure GRAPH_USERNAME XXXXX
python -m knitwork configure GRAPH_PASSWORD XXXXX
python -m knitwork configure GRAPH_LOCATION XXXXX
```

## Running the Fragment process

To run the Fragment process which looks for subnodes and synthons for a given set of fragments/molecules in an SDF:

```
python -m knitwork fragment INPUT_SDF
```
