function [node] = TIDANSE_tree(node,node_update)
% The function performs the T-DANSE algorithm with internal compression,
% essentially the same as DANSE_sum_ic
%
% Syntax:  [node] = t_DANSE_seq(node,node_update)
%
% Inputs:
%   node            -   node structure from network_gen_tree.m
%   node_update     -   which node is updating during this iteration
%
% Outputs:
%   node            -   node structure that contains the cost at each
%                       iteration of the T-DANSE algorithm as well as the
%                       optimal filters
%
% Example:
%    [node] = t_DANSE_seq(node,n)
%
% Other m-files required: none
% Subfunctions: node
% MAT-files required: none
%
% See also: network_gen_tree

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Aug. 2014; Last  revision: 13-Oct-2014

nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;
[node.cost] = deal(0);

for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).P'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).P'*node(ii).ss_noise')';
end

node = TIDANSE_rooted_ff(node,node_update);
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
node = TIDANSE_filt_update(node,node_update);
%%  Calculate cost at each node
for ii=1:nb_nodes
    z_x_seq = node(ii).zx_sum - node(ii).loc_zx;
    z_n_seq = node(ii).zn_sum - node(ii).loc_zn;
    % cost at node during current iteration
    temp_filt = [node(ii).loc_filt_coeff' node(ii).gkq(1).coeff'];
    
    % cost at node during current iteration
    node(ii).cost(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end

end


