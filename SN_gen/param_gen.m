function [sim_param,DANSE_param] = param_gen()
% number of desired sources (also dimension of DANSE)
sim_param.nb_desired_src = 2;
% number of correlated noise sources
sim_param.nb_noise_src = 4;   
% number of nodes
sim_param.nb_nodes = 15;
% length of signal
sim_param.sig_len = 10000;
% size of environment (x-y plane)
sim_param.env_size = [5 5];                   


% Dimension of DANSE  equal to number of desired sources
DANSE_param.dim = sim_param.nb_desired_src;
% max iterations before stopping DANSE algorithms (note algoritm may not
% have converged when max iter has been reached)
DANSE_param.max_iter = 2500;  
% threshold for when to stop algorithms, i.e., when convergence is met
DANSE_param.thresh = 1e-5;  