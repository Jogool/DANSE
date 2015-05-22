function [node,source,noise,white_noise_var] = ...
    network_gen(DANSE_param,sim_param)
% generate random WSN
%
% The function generates a set of nodes randomly placed in a sensing
% environment. 
%
% Syntax:  [node,source,noise,white_noise_var] = network_gen()
%
% Inputs: 
%   DANSE_param     -   DANSE parameters
%   sim_param       -   simulation settings
%
% Outputs:
%   node            -   contains the generated nodes in a structure
%   source          -   raw source signals
%   noise           -   raw noise signals
%   white_noise_var -   uncorrelated additive noise
%
% Example:
%    [node,source,noise,white_noise_var] = network_gen(DANSE_param,sim_param)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: batch_run

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 22-May-2015

% hardcoded parameters
node(sim_param.nb_nodes) = struct;
% generate sources at random locations 
for ii = 1:sim_param.nb_desired_src
    source(ii).nb = ii;                                 % source number
    source(ii).pos = sim_param.env_size.*[rand([1,2])]; % source position
    source(ii).signal = -0.5+rand(sim_param.sig_len,1); % source signal   
end

% generate noise sources at random locations
for ii = 1:sim_param.nb_noise_src
    noise(ii).nb = ii;                                  % noise source number
    noise(ii).pos = sim_param.env_size.*[rand([1,2])];  % noise source position
    noise(ii).signal = -0.5+rand(sim_param.sig_len,1);  % noise signal
end

white_noise_var = mean(var(cat(2,source.signal)))/2;    % representative of sensor noise (uncorrelated between nodes)

% generate nodes at random locations
for ii = 1:sim_param.nb_nodes
    node(ii).nb = ii;                                   % node number
    node(ii).sensors = 5;                               % number of sensors on node
    node(ii).ss_clean = zeros(sim_param.sig_len,node(ii).sensors); % pre-allocate clean source signals
    node(ii).ss_noise = zeros(sim_param.sig_len,node(ii).sensors); % pre-allocate additive noise signals
    node(ii).pos = sim_param.env_size.*[rand([1,2])];              % node position
    
    % this is the sensor array for the node (can customize sensor positions) 
    r_sen = 0.1;                                    % sensor distance from center of node (equispaced) now 10 centimeters
    circ_pos = linspace(0,2*pi,node(ii).sensors+1); % position of sensors around node

    for jj = 1:node(ii).sensors
        node(ii).sensor(jj).pos(1) = node(ii).pos(1)+r_sen*cos(circ_pos(jj));   % x-axis sensor position on node
        node(ii).sensor(jj).pos(2) = node(ii).pos(2)+r_sen*sin(circ_pos(jj));   % y-axis sensor position on node
        for kk = 1: size(source,2)
            d = norm(source(kk).pos - node(ii).sensor(jj).pos);                % calculate attenuation factor based on Euclidean distance
            node(ii).steering(jj,kk) = d;                                      % steering matrix (can be used to check G-coefficients in DANSE algorithm)
            node(ii).ss_clean(:,jj) = node(ii).ss_clean(:,jj)+d*source(kk).signal;
        end
        for kk = 1:size(noise,2)
            d = norm(noise(kk).pos - node(ii).sensor(jj).pos);
            node(ii).ss_noise(:,jj) = node(ii).ss_noise(:,jj)+d*noise(kk).signal+sqrt(white_noise_var)*randn(sim_param.sig_len,1);
        end
    end
    
    % generate local filter coeff and Gkqs for all other nodes
    node(ii).loc_filt_coeff = -1+2*rand(node(ii).sensors,DANSE_param.dim);
    idx = find(ii ~= 1:sim_param.nb_nodes);
    for jj = idx
        node(ii).gkq(jj).coeff = zeros(DANSE_param.dim,DANSE_param.dim);
    end 
    % generate broadcast filter for TI-DANSE
    node(ii).P = eye(size(node(ii).loc_filt_coeff));
end

