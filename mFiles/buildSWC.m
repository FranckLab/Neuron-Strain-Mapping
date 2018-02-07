function swc = buildSWC(tree)
%
% buildSWC(tree)
%
% Create an swc array from a tree structure
%
% INPUTS
%------------------
% tree:          Tree structure stored in a structure class. 
%--------            
% OUTPUTS            
% -----------------
% swc:              Array containing the following columns: [n T x y z R
%                   If multiple model neurons exist, swc will be a cell
%                   array
%
% n: integer label that identifies the current point and increments by one 
%    from one line to the next.
% T: integer representing the type of neuronal segment, such as soma, axon, 
%    apical dendrite, etc.
%    0 = undefined
%    1 = soma (root)
%    2 = axon
%    3 = dendrite
%    4 = apical dendrite
%    5 = fork point
%    6 = end point
%    7 = custom
% x, y, z: cartesian coordinates of each node.
% R: radius at that node.
% P: parent (the integer label) of the current point or -1 to indicate an 
%    origin (soma).

% Place tree structure information into corresponding columns
nN = length(tree.X);
swc = zeros(nN,7);
swc(:,1)    = 1:nN;
swc(:,3)    = tree.X;
swc(:,4)    = tree.Y;
swc(:,5)    = tree.Z;
swc(:,6)    = tree.D / 2; % radius

% Find the parent/child of each node. See Trees Toolbox for documentation
% about tree.dA
PIdx = find(tree.dA);
[parent, child] = ind2sub(size(tree.dA),PIdx);
P = zeros(nN,1);
P(parent) = child;
P(1) = -1;
swc(:,7) = P;

% Set type based on number of connections
%    0 = undefined
%    1 = soma (root)
%    2 = axon
%    3 = dendrite
%    4 = apical dendrite
%    5 = fork point
%    6 = end point
%    7 = custom
swc(1,2) = 1;
for i = 2:length(child)
    conn = sum(i == child);
    switch conn
        case 0, swc(i,2) = 6; % end point
        case 1, swc(i,2) = 3; % dendrite
        otherwise,
            swc(i,2) = 5;
    end
end
end