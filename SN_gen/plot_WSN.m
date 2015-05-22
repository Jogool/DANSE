function [] = plot_WSN(node,A,source,noise)
% plot location of nodes, speech and noise sources
%
% Syntax:  [] = path_find(node,A,source,noise)
%
% Inputs:
%   node            -   contains node information in structure format
%   A               -   adjacency matrix of tree
%   source          -   contains source information in structure format
%   noise           -   contains noise information in structure format
%
% Outputs: none              
% Example:
%    [] = plot_WSN(node,A,source,noise)
%
% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% December 2014; Last revision: 22-May-2015
clf
axis([0 5 0 5]);
hold on

for ii = 1:size(source,2)
    plot(source(ii).pos(1),source(ii).pos(2),'s','MarkerFaceColor','red', ...
        'MarkerEdgeColor','red')
end

for ii = 1:size(noise,2)
    plot(noise(ii).pos(1),noise(ii).pos(2),'d','MarkerFaceColor','blue', ...
        'MarkerEdgeColor','blue')
end

for ii = 1:size(node,2)
    plot(node(ii).pos(1),node(ii).pos(2),'o','MarkerFaceColor','black', ...
        'MarkerEdgeColor','black')
    text(node(ii).pos(1)+0.1,node(ii).pos(2)+0.1,num2str(ii))
end

% plot ad-hoc connections
for ii = 1:size(node,2)
    nb_conn = find(A(ii,:));
    for jj = nb_conn
        h = line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)]);
        set(h,'Color','red')
        set(h,'LineWidth',1)
        set(h,'LineStyle','--')
    end

% plot tree connections
    nb_conn = node(ii).tree_conn';
    for jj = nb_conn
        h = line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)]);
        set(h,'Color','black')
        set(h,'LineWidth',2)
    end 
end
drawnow
