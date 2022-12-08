Noting down what is currently modified by hand in a kinetic sweep:

Always:
In Snakefile:
base dir on line 4
base dir in rule copy_baseline
In scripts/postproc.py:
Copy resources/grids.json from corresponding seeds
(n.b., should probably generate named dirs from seeds and templates)
sites list in postproc, see #TODO

Potentially:
kinetic base and sweeps
timestep send to postproc in classify colonies
cell size hardcoded in classify colonies

Probably only per machine/rarely:
nufeb location

