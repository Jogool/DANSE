function [] = batch_run()
% Syntax:  [cout] = batch_run()
%
% Inputs:  none        - for options see param_gen
%                                                         
% Outputs: none  

% Other m-files required: centralized (see also required folders)
% Folders required: DANSE,TDANSE,TIDANSE,SN_gen)
%
% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% Aug. 2014; Last revision: 22-May-2015

%% hardcoded parameters
[sim_param,DANSE_param] = param_gen;   

plot_on = 1;        % 1(0) - show (do not show) network

% generate random network and TDANSE updating order 
if plot_on
    close all
    [node,source,noise,wnv] = network_gen(DANSE_param,sim_param);
    [node,A,updateorder] = construct_tree(node);
    plot_WSN(node,tril(A),source,noise)
else
    [node,~,~,~] = network_gen(DANSE_param);
    [node,~,updateorder] = construct_tree(node);
end

% find centralized solution
disp('Centralized')
[node] = centralized(node,sim_param.nb_desired_src);

% store original coefficients, this are loaded before every instance of a
% DANSE algorithm, so that the local filters all start at the same value
org_node = node;

% find inital cost (before any updates)
node=init_cost_fc(node,sim_param.nb_nodes,DANSE_param.dim);
cost = sum(cat(1,node.cost));
node = org_node;
%% DANSE
fprintf('\n')
disp('DANSE')
reverseStr = '';

node_update = updateorder(1);
cost_sum_DANSE = cost;
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,DANSE_param.thresh),ge(ii,DANSE_param.max_iter))
    [node] = DANSE(node,node_update,sim_param.nb_nodes,DANSE_param.dim);
    cost_sum_DANSE = [cost_sum_DANSE sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');

    ii = ii + 1;  
    node_update=rem(node_update,sim_param.nb_nodes)+1;    
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
node = org_node;

% find inital cost (before any updates)
node=init_cost_t(node,updateorder(1),sim_param.nb_nodes,DANSE_param.dim);
cost = sum(cat(1,node.cost));
node = org_node;

%% T-DANSE
fprintf('\n')
reverseStr = '';
disp('TDANSE')
node_update = updateorder(1);
cost_sum_TDANSE = cost;
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,DANSE_param.thresh),ge(ii,DANSE_param.max_iter))
    [node] = TDANSE(node,node_update,sim_param.nb_nodes,DANSE_param.dim);
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

% find inital cost (before any updates)
node=init_cost_ti(node,sim_param.nb_nodes,DANSE_param.dim);
cost = sum(cat(1,node.cost));
cost_sum_TIDANSE_fc = cost;
node = org_node;
%% TI-DANSE
fprintf('\n')
reverseStr = '';
disp('TI-DANSE FC')
for ii = 1:sim_param.nb_nodes
    [node(ii).gkq] = deal([]);
    node(ii).gkq(1).coeff = zeros(DANSE_param.dim,DANSE_param.dim);
end
node_update = updateorder(1);
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,DANSE_param.thresh),ge(ii,DANSE_param.max_iter))
    [node] = TIDANSE_fc(node,node_update,sim_param.nb_nodes,DANSE_param.dim);
    cost_sum_TIDANSE_fc = [cost_sum_TIDANSE_fc sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');

    ii = ii + 1;
    node_update=rem(node_update,sim_param.nb_nodes)+1;
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
node = org_node;

%% TI-DANSE tree
fprintf('\n')
reverseStr = '';
disp('TI-DANSE T')
for ii = 1:sim_param.nb_nodes
    [node(ii).gkq] = deal([]);
    node(ii).gkq(1).coeff = zeros(DANSE_param.dim,DANSE_param.dim);
end

node_update = updateorder(1);
cost_sum_TIDANSE_tree = cost;
ii = 1;
tot_diff = inf;
% check if we have met either condition
while ~or(lt(tot_diff,DANSE_param.thresh),ge(ii,DANSE_param.max_iter))
    [node] = TIDANSE_tree(node,node_update,sim_param.nb_nodes,DANSE_param.dim);
    cost_sum_TIDANSE_tree = [cost_sum_TIDANSE_tree sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) -  ...
        cellfun(@(x) x(end), {node.cost})');
    
    ii = ii + 1;
    node_update=rem(node_update,sim_param.nb_nodes)+1;
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n')
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
box on
end




