function [node] = TIDANSE_fc(node,node_update)
[node.cost] = deal(0);
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;

for ii = 1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    node(ii).loc_zx = (node(ii).P'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).P'*node(ii).ss_noise')';
end

%% Contstruct summed signal
z_x_seq = [node.loc_zx];
z_n_seq = [node.loc_zn];
for jj = 1:dim_DANSE;
    zx_sum(:,jj) = sum(z_x_seq(:,jj:dim_DANSE:end),2);
    zn_sum(:,jj) = sum(z_n_seq(:,jj:dim_DANSE:end),2);
end

for ii = 1:nb_nodes
    node(ii).zx_sum = zx_sum;
    node(ii).zn_sum = zn_sum;
end
% update external filter coefficients at root node
node = TIDANSE_filt_update(node,node_update);

%%  Calculate cost at each node
for ii=1:nb_nodes
    z_x_seq = node(ii).zx_sum - node(ii).loc_zx;
    z_n_seq = node(ii).zn_sum - node(ii).loc_zn;
    
    temp_filt = [node(ii).loc_filt_coeff' node(ii).gkq(1).coeff'];
    
    % cost at node during current iteration
    node(ii).cost(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end
