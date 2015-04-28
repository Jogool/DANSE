function [total_conv] = batch_run()
% Syntax:  [cout] = batch_run()
%
% Inputs:   none - for options change hardcoded parameters
%                                                         
%
% Outputs: total_conv       - store the summed LS cost of all nodes
%                             for all algorithms

% Other m-files required: network_gen_tree, centralized, folders
% (DANSE,TDANSE,TIDANSE_fc,TIDANSE_tree)
% MAT-files required: none
%
% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% Aug. 2014; Last revision: 23-Apr-2014

%% hardcoded parameters
% number of desired sources (also dimension of DANSE)
DANSE_param.desired_sources = 2;  

% number of nodes
DANSE_param.nb_nodes = 10;       

% number of sensors per node (assumed same across all nodes compression 
%ratio of (DANSE_param.sensors+1)/DANSE_param.desired_sources)
DANSE_param.sensors = DANSE_param.desired_sources + 1; 

% number of correlated noise sources
DANSE_param.noise_sources = 4;              


plot_on = 1;        % 1(0) - show (do not show) network

% max iterations before stopping DANSE algorithms 
%note algoritm may not have converged when max iter has been reached
max_iter = 1000;  

 % threshold for when to stop algorithms, i.e., when convergence is met
thresh = 1e-5;     

% output 
total_conv = zeros(5,max_iter); % see header for description

% generate random network and TDANSE updating order 
if plot_on
    close all
    [node,source,noise,wnv] = network_gen(DANSE_param);
    [node,A,updateorder] = construct_tree(node);
    plot_WSN(node,A,source,noise)
else
    [node,~,~,~] = network_gen(DANSE_param);
    [node,~,updateorder] = construct_tree(node);
end

% find centralized solution
disp('Centralized')
[node] = centralized(node);

% store original coefficients, this is loaded before every instance of a
% DANSE algorithm, so that the local filters all start at the same value
org_node = node;

% DANSE
fprintf('\n')
disp('DANSE')
reverseStr = '';

node_update = updateorder(1);
cost_sum_DANSE = [];
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,thresh),ge(ii,max_iter));
    [node] = DANSE(node,node_update);
    cost_sum_DANSE = [cost_sum_DANSE sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');

    ii = ii + 1;  
    node_update=rem(node_update,DANSE_param.nb_nodes)+1;    
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

node = org_node;
% T-DANSE
fprintf('\n')
reverseStr = '';
disp('TDANSE')
node_update = updateorder(1);
cost_sum_TDANSE = [];
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,thresh),ge(ii,max_iter));
    [node] = TDANSE(node,node_update);
    cost_sum_TDANSE = [cost_sum_TDANSE sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');
    
    ii = ii + 1;
    node_update=updateorder(rem(ii,numel(updateorder))+1);
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

node = org_node;
% TI-DANSE
fprintf('\n')
reverseStr = '';
disp('TI-DANSE FC')
node=org_node;
for ii = 1:DANSE_param.nb_nodes
    temp = [node(ii).gkq.coeff];
    [node(ii).gkq] = deal([]);
    node(ii).gkq(1).coeff = temp(1:DANSE_param.desired_sources,1:DANSE_param.desired_sources);
    node(ii).P =  node(ii).loc_filt_coeff / node(ii).gkq(1).coeff;
end
cost_sum_TIDANSE_fc = [];
node_update = updateorder(1);
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,thresh),ge(ii,max_iter));
    [node] = TIDANSE_fc(node,node_update);
    cost_sum_TIDANSE_fc = [cost_sum_TIDANSE_fc sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');

    ii = ii + 1;
    node_update=rem(node_update,DANSE_param.nb_nodes)+1;
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

node=org_node;
% TI-DANSE tree
fprintf('\n')
reverseStr = '';
disp('TI-DANSE T')
for ii = 1:DANSE_param.nb_nodes
    temp = [node(ii).gkq.coeff];
    [node(ii).gkq] = deal([]);
    node(ii).gkq(1).coeff = temp(1:DANSE_param.desired_sources,1:DANSE_param.desired_sources);
    node(ii).P =  node(ii).loc_filt_coeff / node(ii).gkq(1).coeff;
end
node_update = updateorder(1);
cost_sum_TIDANSE_tree = [];
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,thresh),ge(ii,max_iter));
    [node] = TIDANSE_tree(node,node_update);
    cost_sum_TIDANSE_tree = [cost_sum_TIDANSE_tree sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) -  ...
        cellfun(@(x) x(end), {node.cost})');
    

    ii = ii + 1;
    node_update=rem(node_update,DANSE_param.nb_nodes)+1;
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

if plot_on
    figure
    hold on
    loglog(cost_sum_DANSE)
    loglog(cost_sum_TDANSE,'-xm')
    loglog(cost_sum_TIDANSE_fc,'-or')
    loglog(cost_sum_TIDANSE_tree,'--dk')
    axis tight

h =  refline(0,sum([node.cost_cent]));
set(h,'LineStyle','--');

a = get(gca,'YLim');
set(gca,'YLim',[sum([node.cost_cent]) - sum([node.cost_cent])*.1 a(2)])
legend('DANSE', 'T-DANSE', 'TI-DANSE (FC)','TI-DANSE (T)', 'Optimal');
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Iteration')
ylabel('Sum of LS cost for all nodes (dB)')
end

total_conv(1,1:length(cost_sum_DANSE)) = cost_sum_DANSE;
total_conv(2,1:length(cost_sum_TDANSE)) = cost_sum_TDANSE;
total_conv(3,1:length(cost_sum_TIDANSE_fc)) = cost_sum_TIDANSE_fc;
total_conv(4,1:length(cost_sum_TIDANSE_tree)) = cost_sum_TIDANSE_tree;
total_conv(5,1) = sum([node.cost_cent]);
fprintf('\n')



