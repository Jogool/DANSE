DANSE
=====
Code for the DANSE algorithms

batch_run.m is a wrapper script that will run all the DANSE algorithms and compare them.  Note that the parameters are hardcoded, which are easily accessible.

! Before running make sure the subfolders are added to your path.

WSN-generation
==============

These files generate a random WSN and node structure that can be used for the provided DANSE algorithms

* network_gen.m     -  generates the node strucutre 
* construct_tree.m  - generates an ad-hoc network and then prune to tree (note ad-hoc connections are still stored at a node)
* path_find.m       - generates the updating order for the TDANSE algorithm (note that it is only a sufficient not                         necessary condition for the TDANSE algortihm to converge)
* plot_WSN.m        - used for plotting (optional)

DANSE_param.mat contains a sample structure that can be used to generate the network using network_gen.m
node.mat contains a sample node structure

