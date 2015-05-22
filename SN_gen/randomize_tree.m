function [node,A] = randomize_tree(node)
% constrcut a MST based on current node positions
%
% The function generates a MST based on current node positions
%
% Syntax:  [node] = network_gen_tree(node)
%
% Inputs:
%   node            -   contains node information in structure format
%
% Outputs:
%   node            -   contains node with updated tree connections 
%                       in structure format
%
% Example:
%    [node] = construct_tree(node)
%
% Other m-files required: none
% Subfunctions: path_find
% MAT-files required: none
%
% See also: 

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 10-Oct-2014

% increase transmission radius of each node until network is connected
nb_nodes = size(node,2);
A = zeros(nb_nodes);
for ii = 1:nb_nodes
    idx = node(ii).conn;
    A(idx,ii) = 1;
    A(ii,idx) =1 ;
end
A_tril = tril(A);
euc_weight = rand(length(find(A_tril)),1);

[rows,cols] = find(A_tril);


% Bio-graphs
UG = sparse(rows,cols,euc_weight,nb_nodes,nb_nodes);
%view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))

[ST,~] = graphminspantree(UG);
A_mst = full(ST)+full(ST)';    % adjanceny matrix of tree topology

% populate tree connections of each node
for ii = 1:nb_nodes
    node(ii).tree_conn = find(A_mst(:,ii));
end
