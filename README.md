# par-max-dicut

## Overview

This repository contains par-max-dicut: A shared memory framework to compute high quality Max-Dicuts in graphs. This framework uses a graph partitioning approach to partition a graph into k subgraphs and then, on each of the k subgraphs in parallel, runs a sequential approximation algorithm for Max-Dicut. In a final step the local cuts are merged into a local solution. In an optional step the computed cut is further optimized by performing a local search.

## Dependencies

This framework has the following dependencies:
* A compiler supporting C++20
* OpenMP
* Intel Thread Building Blocks library (TBB)
* Optional:
    * Gurobi (To solve Linear Programs
    * Mosek (To solve Semidefinite Programs)

## Building

To build this framework you need to perform the following steps in the root directory:

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
    
## Running benchmarks

The following is an exemplary command to use our framework:

```
./benchmark -r 1 --partitioner part-kaminpar --seqmc seq-goemans --merger merge-mc-tree[seq-goemans,256] -t 8 -K 32 -o --type [mtx|edges|rudy|metis] <file>
```

This command uses the partitioner `part-kaminpar` to partition the graph `<file>` (the graph format is specified by the parameter `--type`) into 32 subgraphs and runs on each subgraph the algorithm `seq-goemans`. The algorithm `merge-mc-tree[seq-goemans,256]` is used to merge the solutions of each subgraph into a global solution by using the algorithm `seq-goemans`. The flag `-o` indicates that the solution is optimized by using a local search in a post-processing step. Overall the framework uses 8 threads and is executed once (specified by the parameter `-r 1`).

The output is in a format that can be processed by Bingmann's [sqlplot-tools](https://github.com/bingmann/sqlplot-tools). An exemplary output can is shown below:

```
RESULT dir=321689 rev=321689 undir=643378 time=13708 algo=par-meta[part-kaminpar,seq-goemans,merge-mc-tree[seq-goemans,256]] K=32 threads=8 run=0 n=11586 m=1136618 file=<file> partitioning=263 seqmc=13169 merge=44 optimize=231
```

Here `dir` refers to the computed directed cut, `rev` to the reversed directed cut and `undir` to the undirected cut. `time` is the total runtime of our framework, `partitioning` the runtime for the partitioning, `seqmc` the runtime to compute all sequential cuts for each subgraph, `merge` the runtime to merge the local solutions and `optimize` is the runtime for the local search. `algo` specifies the configuration of our framework, `K` is the number of subgraphs the graph `file` (with `n` nodes and `m` edges) is partitioned into. `threads` specifies the number of used threads and `run` the number of execution of our framework.
