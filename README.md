# moving-points-polynomial

This is a program to calculate asymptotics of the number of moving points in a strongly connected directed graph.

The repository consists of the following files:

- `hamiltonian.h` -- implementation of class `HamiltonianGraph`, which checks directed graph for being hamiltonian;
- `kosaraju.h` -- implementation of class `SCCGraph`, which checks directed graph for being strongly connected;
- `matrix.h` -- implementation of class `Matrix` for working with matrices over the field of rational numbers;
- `input.txt` -- this is the file where you define your graph; currently contains a bunch of example graphs;
- `main.cpp` -- main file, which computes the power and leading coefficient of the desired asymptotics;

How to run the code:

In the `main` function there are five constant variables which determine the behavior of the program:

- `graph_type` -- determines which type of graph format is being given the `input.txt` file (see below for more);
- `only_print` -- if set `true`, only prints the polynomial and a bunch of statistics of the graph;
- `step_size` -- determines how often to print logs of the program;
- `mem_last` -- smoothing parameter for calculating actual statistics;
- `max_time` -- actual working time of the program.

There are five types of input formats for the graph.

If `graph_type = 1`, then it has the following format:

```cpp
n m
u_1 v_1 w_1
u_2 v_2 w_2
...
u_m v_m w_m
src
```

where `n` is the number of vertices in the graph, `m` is the number of edges in the graph, `u_i v_i w_i` is the `i`-th edge in the graph, which means an edge from the vertex `u_i` to `v_i` of length `w_i`, and `src` is the starting point of the graph dynamics.

If `graph_type = 2`, then the graph has the following format:

```cpp
n m
u_1 v_1
u_2 v_2
...
u_m v_m
src
```

where all of the variables mean the same as in the previous case, but this time the edge lengths are not specified in the input and are randomly generated within the program.

If `graph_type = 3`, then the graph has the following format:

```cpp
n m
u_1 v_1 p_1
u_2 v_2 p_2
...
u_m v_m p_m
src
```

where all of the variables mean the same as in the first case, but this time `p_i` is a prime number and it means that the length of the `i`-th edge is going to be `sqrt(p_i)`.

If `graph_type = 4`, then the graph is not read from the input, but completely generated within the program. In this case, a strongly connected directed (but not hamiltonian) graph is generated.

If `graph_type = 5`, then the graph is not read from the input, but completely generated within the program. In this case, a strongly connected hamiltonian directed graph is generated.

While running, the program periodically outputs three numbers:

```cpp
cur_time cur_points p_exact p_lb p_ub
```

which mean:

- `cur_time` -- time that has passed from the start of the point movement;
- `cur_points` -- exact number of the points in the graph;
- `p_exact` -- exact approximation of the number of points;
- `p_lb` -- lower bound on the number of points;
- `p_ub` -- upper bound on the number of points.
