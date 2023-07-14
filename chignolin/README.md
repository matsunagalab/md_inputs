# chignolin input files for molecular dynamics simulation.

Repository for teaching beginners what is molecular dynamics simulation

Let's find the native structure (1uao.pdb, beta-sheet-like) of chignolin starting from an alpha-helix-like structure. 

## Jupyter notebook

```
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/work matsunagalab/openmm
```

## Execute a script

```
# On Mac
docker run --rm -v "$PWD":/work matsunagalab/openmm python -u run.py

# On Linux
docker run --privileged -u $(id -u):$(id -g) -e CUDA_VISIBLE_DEVICES="0" --rm -v "$PWD":/work matsunagalab/openmm python -u run.py
```

