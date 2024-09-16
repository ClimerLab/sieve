# sieve
A community detection algorithm based on optimizing the Sieve objective function using IBM ILOG CPLEX to solve an Integer Linear Program.

## Setup
Clone repo

Install IBM ILOG CPLEX

Update Makefile with path to CPLEX and CONCERT directories

Compile using the Makefile by navigating to the root directory and entering: make

## Execute Program
This repo contains 3 programs that are to be used together to optimize the Sieve objective function of a graph. The first program, 'SepComp', identifies all of the seperated components in the original network and creates to files for each compnent. Next, 'S_MIP' is executed on each component to find an optimal solution for the component. Finally, 'CombClustOut' combines the results from each component into an optimal solution of the original network.

The user can call each program seperately or they may use the included shell script.

ORIG_GML - Full file path of the .gml or .graph representing the network

OUTPUT_DIR - Directory to save the output and temporary files. You must include the trailing '/'

RUNTAG - Name included in output files for run

CLEANUP_FLAG - Boolean used to control cleanup of temp files

## Output
S_ObjValue.csv - File containint the component number, Sieve value, and run time of each component

<RUNTAG>_S.out - Optimal community assignment for the nodes in the network

<RUNTAG>.stats - Statistics from SepComp