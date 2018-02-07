function [E, L, X, Y, Z] = edgeTree(tree)
%
% edgeTree(tree)
%
% Calculates the lengths between each node 
%
% INPUTS
%------------------
% tree:             Tree structure. Can be either cell array or structure.
%--------          
% OUTPUTS            
% -----------------
% E:                The index to the direct parent node of each individual 
%                   element in the tree 
% L:                Length between each node
% X,Y,Z:            Cartesian coordinates

% Make sure tree is a cell
if ~iscell(tree), tree = {tree}; end

% Pre-allocate
nSoma = length(tree);
E = cell(nSoma,1);
X = cell(nSoma,1);
Y = cell(nSoma,1);
Z = cell(nSoma,1);
L = cell(nSoma,1);

% For number of neurons
for i = 1:nSoma,
    [E{i}, L{i}, X{i}, Y{i}, Z{i}] = edgeTree_fun(tree{i});
end

if nSoma == 1,
    E = E{1};
    X = X{i};
    Y = Y{i};
    Z = Z{i};
    L = L{i};
end

end

function [E, L, X, Y, Z] = edgeTree_fun(tree)

% Returns the index to the direct parent node of each individual element in
% the tree 
E(:,2) = idpar_tree(tree);

% Simply numbers the nodes
E(:,1) = 1:length(E(:,1));

% Coordinates of each node and parent node
X = [tree.X(E(:,1)), tree.X(E(:,2))];
Y = [tree.Y(E(:,1)), tree.Y(E(:,2))];
Z = [tree.Z(E(:,1)), tree.Z(E(:,2))];

% Difference
dX = X(:,2) - X(:,1);
dY = Y(:,2) - Y(:,1);
dZ = Z(:,2) - Z(:,1);

% Edge length
L = sqrt(dX.^2 + dY.^2 + dZ.^2);

end