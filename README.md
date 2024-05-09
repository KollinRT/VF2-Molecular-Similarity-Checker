# Overview of project
## Background
This project stemmed from exploring algorithms that are utilized for solving Graph Isomorphism problems. I chose to focus on the Graph-Subgraph Isomorphism problem and wanted to explore how to use graph algorithms to determine the largest common substructures between graphs.

## Motivation
Molecules are ubiquitous. They make up drugs, materials, and many more. 

Primary motivation is due to a medicinal chemistry course that was taken and understanding if there was a graph-theoretic approach to determining pharmacophores for pairs or larger of molecules.

## Background on Algorithm
The VF2 algorithm comes from the paper "An Improved Algorithm for Matching Large Graphs" written by L. P. Cordella, P. Foggia, C. Sansone, and M. Vento. The idea behind the algorithm is shown in the following figure taken from the paper.

![VF2Algo](C:\Users\kolli\Desktop\VF2Algo\VF2-Molecular-Similarity-Checker\images\VF2Algo.png)

Which shows us that this Match algorithm scans for common matches and grows the subgraph and if not include then includes backtracking  to help determine the best most connected graph.

## Implementation

The implementation was carried out in Python referencing some packages such as RDKit.

## TODO

\- [ ] Clean up the Jupyter Notebook for examples of utilization beyond the main script
