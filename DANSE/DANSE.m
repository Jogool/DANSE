function [node] = DANSE(node,node_update,nb_nodes,dim_DANSE)
% perform a single iteration of the DANSE algorithm 
%
% Syntax:  [node] = DANSE(node,node_update,nb_nodes,dim_DANSE)
%
% Inputs:
%   node            -   contains node information in structure format
%   node_update     -   which node to update
%   nb_nodes        -   number of nodes
%   dim_DANSE       -   number of broadcast signals per node
%
% Outputs:
%   node            -   contains node with updated filter and cost
%                   
% Example:
%     [node] = DANSE(node,node_update,nb_nodes,dim_DANSE)
%
% Other m-files required: DANSE_filt_update

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 22-May-2015
[node.cost_DANSE] = deal(0);
% generate broadcast signals
for ii = 1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end
% update external filter coefficients at root node
node = DANSE_filt_update(node,node_update,nb_nodes,dim_DANSE);

%%  Calculate cost at each node
for ii=1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    z_x_seq = [node(idx).loc_zx];
    z_n_seq = [node(idx).loc_zn];
    
    gkq_coeff = cat(1,node(ii).gkq.coeff);
    temp_filt = [node(ii).loc_filt_coeff' gkq_coeff'];
    
    % cost at node during current iteration
    node(ii).cost = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end