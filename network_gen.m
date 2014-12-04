function [node,source,noise,white_noise_var] = network_gen(DANSE_param)
% generate random WSN
%
% The function generates a set of nodes randomly placed in a sensing
% environment. 
%
% Syntax:  [node,source,noise,white_noise_var] = network_gen()
%
% Inputs: none (passed via global parameters)
%
% Outputs:
%   node            -   contains the generated nodes in a structred format
%   source          -   raw source signals
%   noise           -   raw noise signals
%   white_noise_var -   uncorrelated additive noise
%
% Example:
%    [node,source,noise,white_noise_var] = network_gen()
%
% Other m-files required: batch_run
% Subfunctions: none
% MAT-files required: none
%
% See also: batch_run

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 13-Oct-2014


% hardcoded parameters
length_of_signal = 10000;       % length or source signals
r_sen = 0.1;                    % sensor distance from center of node (equispaced) now 10 centimeters
circ_pos = linspace(0,2*pi,DANSE_param.sensors+1); % position of sensors around node
node(DANSE_param.nb_nodes) = struct;
[node.dimDANSE] = deal(DANSE_param.desired_sources);
% generate sources at random locations 
for ii = 1:DANSE_param.desired_sources
    source(ii).nb = ii;                     % source number
    source(ii).pos  = 5*rand(2,1);          % source position (assuming 5x5 room)
    source(ii).signal = -0.5+rand(length_of_signal,1);  % source signal   
end

% generate noise sources at random locations
for ii = 1:DANSE_param.noise_sources
    noise(ii).nb = ii;                      % noise source number
    noise(ii).pos  = 5*rand(2,1);           % noise source position (assuming 5x5 room)
    noise(ii).signal = -0.5+rand(length_of_signal,1);   % noise signal
end

white_noise_var = mean(var(cat(2,source.signal)))/2;    % representative of sensor noise (uncorrelated between nodes)

% generate nodes at random locations
for ii = 1:DANSE_param.nb_nodes
    node(ii).nb = ii;                                   % node number
    node(ii).sensors = DANSE_param.sensors;             % number of sensors on node
    node(ii).ss_clean = zeros(length_of_signal,node(ii).sensors); % pre-allocate clean source signals
    node(ii).ss_noise = zeros(length_of_signal,node(ii).sensors); % pre-allocate additive noise signals
    node(ii).pos = 5*rand(2,1);                         % node position (assuming 5x5 room)
    
    for jj = 1:node(ii).sensors
        node(ii).sensor(jj).pos(1) = node(ii).pos(1)+r_sen*cos(circ_pos(jj));   % x-axis sensor position on node
        node(ii).sensor(jj).pos(2) = node(ii).pos(2)+r_sen*sin(circ_pos(jj));   % y-axis sensor position on node
        for kk = 1: size(source,2)
            d = norm(source(kk).pos - node(ii).sensor(jj).pos');                % calculate attenuation factor based on Euclidean distance
            node(ii).steering(jj,kk) = d;                                       % steering matrix (can be used to check G-coefficients in DANSE algorithm)
            node(ii).ss_clean(:,jj) = node(ii).ss_clean(:,jj)+d*source(kk).signal;
        end
        for kk = 1:size(noise,2)
            d = norm(noise(kk).pos - node(ii).sensor(jj).pos');
            node(ii).ss_noise(:,jj) = node(ii).ss_noise(:,jj)+d*noise(kk).signal+sqrt(white_noise_var)*randn(length_of_signal,1);
        end
    end
    
    % generate local filter coeff and Gkqs for all other nodes
    node(ii).loc_filt_coeff = -1+2*rand(node(ii).sensors,DANSE_param.desired_sources);
    idx = find(ii ~= 1:DANSE_param.nb_nodes);
    for jj = idx
        node(ii).gkq(jj).coeff = -1+2*rand(DANSE_param.desired_sources,DANSE_param.desired_sources);
    end 
end

