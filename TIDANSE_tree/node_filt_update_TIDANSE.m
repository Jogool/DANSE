%% update filter coefficents of node
function [node] = node_filt_update_TIDANSE(node,ii)
global dim_DANSE
z_x_seq = node(ii).zx_sum - node(ii).loc_zx;
z_n_seq = node(ii).zn_sum - node(ii).loc_zn;


Rxx = [node(ii).ss_clean z_x_seq]'*[node(ii).ss_clean z_x_seq];
Rnn = [node(ii).ss_noise z_n_seq]'*[node(ii).ss_noise z_n_seq];

w_temp  = (Rnn+Rxx) \ Rxx(:,1:dim_DANSE);      
node(ii).loc_filt_coeff = w_temp(1:node(ii).sensors,:);                     
node(ii).gkq(1).coeff = w_temp(node(ii).sensors+1:end,:);

node(ii).P = node(ii).loc_filt_coeff /  ...                                 % store P
    node(ii).gkq(1).coeff;
end