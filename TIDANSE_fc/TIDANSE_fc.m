function [node] = TIDANSE_fc(node,node_update)
global nb_nodes dim_DANSE

[node.cost_DANSE] = deal(0);

for ii = 1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    node(ii).loc_zx = (node(ii).P'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).P'*node(ii).ss_noise')';
end

% update external filter coefficients at root node
node = node_filt_update_TIDANSE(node,node_update);

%%  Calculate cost at each node
for ii=1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    z_x_seq = [node(idx).loc_zx];
    z_n_seq = [node(idx).loc_zn];
    
    for jj = 1:dim_DANSE;
        z_x_sum(:,jj) = sum(z_x_seq(:,jj:dim_DANSE:end),2);
        z_n_sum(:,jj) = sum(z_n_seq(:,jj:dim_DANSE:end),2);
    end
    
    temp_filt = [node(ii).loc_filt_coeff' node(ii).gkq(1).coeff'];
    
    % cost at node during current iteration
    node(ii).cost(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_sum]+[node(ii).ss_noise z_n_sum])')^2;
end
