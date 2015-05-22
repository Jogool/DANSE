function [node] = TIDANSE_tree(node,node_update,nb_nodes,dim_DANSE)
% perform a single iteration of the TI-DANSE algorithm in a tree topology
%
% Syntax:  [node] = TIDANSE_tree(node,node_update,nb_nodes,dim_DANSE)
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
%     [node] = TIDANSE_tree(node,node_update,nb_nodes,dim_DANSE)
%
% Other m-files required: TIDANSE_filt_update

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 22-May-2015
[node.cost] = deal(0);

for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).P'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).P'*node(ii).ss_noise')';
end

node = TIDANSE_rooted_ff(node,node_update,nb_nodes,dim_DANSE);

%% Contstruct summed signal at root and flood signal through WSN
idx = node(node_update).ff_rec;
z_x_seq = [node(idx).ff_zx];
z_n_seq = [node(idx).ff_zn];
for jj = 1:dim_DANSE;
    z_x_sum(:,jj) = sum(z_x_seq(:,jj:dim_DANSE:end),2);
    z_n_sum(:,jj) = sum(z_n_seq(:,jj:dim_DANSE:end),2);
end
zx_sum = node(node_update).loc_zx + z_x_sum;
zn_sum = node(node_update).loc_zn + z_n_sum;
for ii = 1:nb_nodes
    node(ii).zx_sum = zx_sum;
    node(ii).zn_sum = zn_sum;
end

% update node-specific filter coefficients at updated node
node = TIDANSE_filt_update(node,node_update,dim_DANSE);
%%  Calculate cost at each node
for ii=1:nb_nodes
    z_x_seq = node(ii).zx_sum - node(ii).loc_zx;
    z_n_seq = node(ii).zn_sum - node(ii).loc_zn;
    % cost at node during current iteration
    temp_filt = [node(ii).loc_filt_coeff' node(ii).gkq(1).coeff'];
    
    % cost at node during current iteration
    node(ii).cost = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end

end


