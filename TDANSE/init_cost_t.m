function node = init_cost_t(node,root,nb_nodes,dim_DANSE)
for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end

node = TDANSE_rooted_ff(node,root,nb_nodes,dim_DANSE);
node = TDANSE_rooted_df(node,root,nb_nodes,dim_DANSE);

for ii=1:nb_nodes
    idx = node(ii).ff_rec;
    z_x_seq = [node(node(ii).ff_trans).df(ii).zx node(idx).ff_zx];
    z_n_seq = [node(node(ii).ff_trans).df(ii).zn node(idx).ff_zn];
    
    gkq_coeff = [node(ii).gkq(node(ii).ff_trans).coeff node(ii).gkq(idx).coeff];
    gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
    gkq_coeff = cat(1,gkq_coeff{:});
    
    % cost at node during current iteration
    node(ii).cost = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        [node(ii).loc_filt_coeff' gkq_coeff']*...
        ([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end

