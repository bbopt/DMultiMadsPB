# DMulti-MADS (PB, TEB and Penalty variants)

This repertory contains the source code of the DMulti-MADS algorithm for constrained blackbox optimization.
In terms of performance, it is more efficient than the old implementation found in [DMulti-MADS][https://github.com/bbopt/DMultiMadsEB].

> Jean Bigeon, Sébastien Le Digabel and Ludovic Salomon, *Handling of constraints in multiobjective blackbox optimization*

**Warning** : This code has no vocation to be used in industry, see [Nomad](https://www.gerad.ca/nomad) for a more robust implementation of state-of-the-art blackbox method.
It aims at guaranteeing the reproducibility of the experiments described in this work.

## Use

To use DMulti-MADS, Julia >= 1.6 is required. One can test it by typing the following command at the root of the directory.
````
julia> ]

(@v1.6) pkg> activate .

(DMultiMadsPB) pkg> test
````
All tests should pass.

A simple example can be found in the *examples/* folder. One can look also at *./test/madsmodel.jl* for more examples.

## Problems

This folder contains an implementation of all multiobjective benchmark problems used in this article for the Nomad (BiMADS) software.

The algorithm [DFMO][http://www.iasi.cnr.it/~liuzzi/DFL/] is provided with all benchmarks coded in Fortran by the authors. For this reason, it is not given here.
> **Warning** The generation of analytical benchmarks takes a lot of time (around three days and requires more than 40 G of memory on hardware;

For BiMADS, all executables are given with models and nelder-mead search _deactivated_. Uncomment the lines in the *main* function if you need them.

To obtain the real blackbox optimization applications, one can get them at:
- [STYRENE][https://github.com/bbopt/styrene]
- [SOLAR][https://github.com/bbopt/solar]
You can compile them and launch them on the different algorithms with the scripts provided (minus paths to adapt to your machines):
> **Warning** Solving STYRENE and SOLAR for a given solver takes one day.
