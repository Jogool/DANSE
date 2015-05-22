function node = init_cost_ti(node,nb_nodes,dim_DANSE)
%%  Calculate inital broadcast signals
for ii = 1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    node(ii).loc_zx = (node(ii).P'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).P'*node(ii).ss_noise')';
end

%%  Calculate sum of signals
z_x_seq = [node.loc_zx];
z_n_seq = [node.loc_zn];

for jj = 1:dim_DANSE;
    zx_sum(:,jj) = sum(z_x_seq(:,jj:dim_DANSE:end),2);
    zn_sum(:,jj) = sum(z_n_seq(:,jj:dim_DANSE:end),2);
end

%%  Calculate cost at each node
g_temp = zeros(dim_DANSE,dim_DANSE);
for ii=1:nb_nodes
    z_x_seq = zx_sum - node(ii).loc_zx;
    z_n_seq = zn_sum - node(ii).loc_zn;
    
    temp_filt = [node(ii).loc_filt_coeff' g_temp];
    
    % cost at node during current iteration
    node(ii).cost = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end
