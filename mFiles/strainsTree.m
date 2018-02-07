function [E_local, E_mean, K_mean, Ef, tree, Ef_eig] = strainsTree(swc,Einf,location)
%
% strainsTree(swc,Einf,location)
%
% Determine the local and global strains of all model neurons included in
% the swc file
%
% INPUTS
%------------------
% swc:          swc file from Neuron Studio
% Einf:         Farfield strain
% location:     Vertex 'v' (node) or edge 'e' (along the segment)
%--------          
% OUTPUTS            
% -----------------
% If swc file contains more than one model neuron, outputs will be stored
% as cell arrays, one cell for each neuron.
%
% E_local:      Local strain values. If swc file contains more than one
%                   model neuron, E_local will be a cell array
%               E_local(:,1) and E_mean(1) -----> radial strain
%               E_local(:,2) and E_mean(2) -----> radial shear strain
%               E_local(:,3) and E_mean(3) -----> axial strain
%               E_local(:,4) and E_mean(4) -----> axial shear strain
%
% E_mean:       Global strain values stored in the same manner as E_local
% 
% K_mean:       Mean curvature of the neuron
% 
% Ef:           3 x 3 x n array, where n = number of nodes. Each 3 x 3 is 
%                   the rotated strain field at each node.
% 
% tree:         Contains the tree structures of all model neurons
%
% Ef_eig:       Local strain values calculated at the edge (along segment)
%--------          
% NOTICE            
% -----------------
% Reason for NaN values is that endpoints and branch points are not
% considered and set to NaN. This can be changed by using splineTreeEdited
% instead of splineTree (first line of code in local function
% strainTree_local)

if nargin < 3, location = 'v'; end % vertex or edge

% Read in tree structure
[~, tree] = readSWC(swc);

% Make nodes equidistant
tree = resampleTree(tree);

% Extract segments from tree structure
r = findSeg(tree);

% Taubin smoothing 
tree = taubinSmoothing(r,tree);

% Initialize variables
nSoma = length(tree);
E_local = cell(nSoma,1);
Ef = cell(nSoma,1);
E_mean = cell(nSoma,1);
K_mean = cell(nSoma,1);
Ef_eig = cell(nSoma,1);

% Compute local and global strains for every neuron
for i = 1:nSoma
    
    % Reformat tree structure into an array
    swc = buildSWC(tree{i});
    % Compute local metrics
    [E_local{i}, Ef{i}] = strainTree_local(swc, Einf, location);
    % Compute global metrics
    [E_mean{i}, K_mean{i}, Ef_eig{i}] = strainTree_global(tree{i},swc, Einf);
end

% Parse outputs
if nSoma == 1 
    E_local = E_local{1};
    Ef = Ef{1}; 
    E_mean = E_mean{1};
    K_mean = K_mean{1};
    Ef_eig = Ef_eig{1};
    tree = tree{1};
end

end

% Local strains
function [Ef_eig,Ef,t,n,b,K,T] = strainTree_local(swc, Einf, location)

% Frenet frame
[t,n,b,K,T] = splineTree(swc,location);
    
% Align tangent with loading axis e3
f1 = b;
f2 = n;
f3 = t;

% Construct rotation matrix 
R = [reshape(f1(:,1),1,1,[]), reshape(f1(:,2),1,1,[]), reshape(f1(:,3),1,1,[]);
    reshape(f2(:,1),1,1,[]), reshape(f2(:,2),1,1,[]), reshape(f2(:,3),1,1,[]);
    reshape(f3(:,1),1,1,[]), reshape(f3(:,2),1,1,[]), reshape(f3(:,3),1,1,[]);];
Ef = zeros(3,3,size(t,1));

% Rotate far field strain
for i = 1:3
    for j = 1:3
        for k = 1:3
            for m = 1:3
                Ef(i,j,:) = Ef(i,j,:) +  R(i,k,:).*Einf(k,m).*R(j,m,:);
            end
        end
    end
end

% Compute eigenvalues for Ef

% Max radial shear strain
E12 = sqrt(( (Ef(1,1,:) - Ef(2,2,:))/2 ).^2 + Ef(1,2,:).^2);

% Min radial strain
E11min = 0.5*(Ef(1,1,:) + Ef(2,2,:)) - E12;

% Compressive axial strain
E33 = Ef(3,3,:);

% Max axial shear strain
Ei3 = sqrt(Ef(1,3,:).^2 + Ef(2,3,:).^2);

% Reshape
E11min = reshape(E11min,[],1);
E12 = reshape(E12,[],1);
E33 = reshape(E33,[],1);
Ei3 = reshape(Ei3,[],1);

Ef_eig = [E11min E12 E33 Ei3];

end

function [E_mean, K_mean,Ef_eig] = strainTree_global(tree,swc,Einf)

% Use edge location this time 
[Ef_eig, ~, ~, ~, ~, K] = strainTree_local(swc, Einf, 'e');

% Edges and connecting vector
[~, Li] = edgeTree(tree);

% Compute mean strain
E_mean = zeros(1,size(Ef_eig,2));
L = nansum(Li(:));
for i = 1:size(Ef_eig,2)
    E_mean(1,i) = nansum(Ef_eig(:,i).*Li)/L;
end

% Compute mean curvature
K_mean = nansum(K.*Li)/L;

end