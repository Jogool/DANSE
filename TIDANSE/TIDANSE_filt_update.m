function [node] = TIDANSE_filt_update(node,ii,dim_DANSE)

z_x_seq = node(ii).zx_sum - node(ii).loc_zx;
z_n_seq = node(ii).zn_sum - node(ii).loc_zn;

for jj = 1:dim_DANSE;
    z_x_sum(:,jj) = sum(z_x_seq(:,jj:dim_DANSE:end),2);
    z_n_sum(:,jj) = sum(z_n_seq(:,jj:dim_DANSE:end),2);
end

Rxx = [node(ii).ss_clean z_x_sum]'*[node(ii).ss_clean z_x_sum];
Rnn = [node(ii).ss_noise z_n_sum]'*[node(ii).ss_noise z_n_sum];

w_temp  = (Rnn+Rxx) \ Rxx(:,1:dim_DANSE);      
node(ii).loc_filt_coeff = w_temp(1:node(ii).sensors,:);                     
node(ii).gkq(1).coeff = w_temp(node(ii).sensors+1:end,:);

node(ii).P = node(ii).loc_filt_coeff /  ...                                 % store P
    node(ii).gkq(1).coeff;
end