function [swc, tree]= readSWC(filename)
%
% readSWC(filename)
%
% Read in the tree structures from the swc file. 
%
% INPUTS
%------------------
% filename:         swc file from Neuron Studio specified as a string
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
%
% tree:                 structure with the following fields: dA, X, Y, Z, 
%                       D,R,rnames. If multiple model neurons exist, tree
%                       will be a cell array
%
% dA: adjacency matrix. See page 5 of trees toolbox documentation.
% X,Y,Z: cartesian coordinates of each node
% D: diameter at that node
% R: integer representing the type of neuronal segment. 
% rnames: which kind of neuronal segments appear in currect tree structure
% 

% Create swc array
swc = parseInputs(filename);

% Create tree structure
tree = cell(length(swc),1);
for i = 1:length(swc),    
    if isstruct(swc{i})
        tree{i} = swc{i};
        swc{i} = buildSWC(tree{i});
    else
        tree{i} = buildTree(swc{i});
    end   
end

% if only one soma in swc
if iscell(swc) && length(swc) == 1,
    swc = swc{1};
    tree = tree{1};
end

end

% Create swc array
function swc = parseInputs(filename)
    
    % Use v3d library to read in file
    swc0 = load_v3d_neuron_file(filename);
    
    % Number of soma
    idxSoma = find(swc0(:,7) == -1)';
    nSoma = numel(idxSoma);
    
    % Mark somas with -1 
    swc = cell(nSoma,1);
    for i = 1:nSoma,
        if nSoma > 1 && i < nSoma
            swc{i} = swc0(idxSoma(i):idxSoma(i+1)-1,:);
        else
            swc{i} = swc0(idxSoma(i):end,:);
        end
        swc{i}(swc{i}(:,2) == 0,2) = 6;
        
        if nSoma > 1
            swc{i}(2:end,7) = swc{i}(2:end,7) - swc{i}(1,1) + 1;
            swc{i}(:,1) = swc{i}(:,1) - swc{i}(1,1) + 1;
        end
        swc{i} = swc{i}(:,1:7);
        swc{i}(1,7) = -1; 
    end
end

function tree = buildTree(swc)

N = size (swc, 1);

% vector containing index to direct parent
P   = swc (:, 7);

dA = sparse(N, N);
for ward = 2 : N;
    dA(ward, P (ward)) = 1;
end
tree.dA = dA;
tree.X  =            swc (:, 3);     % X-locations of nodes on tree
tree.Y  =            swc (:, 4);     % Y-locations of nodes on tree
tree.Z  =            swc (:, 5);     % Z-locations of nodes on tree
tree.D  =            swc (:, 6) * 2; % local diameter values of nodes on tree
[i1, ~, i3] = unique (swc (:, 2));
tree.R  = i3;
tree.rnames = cellstr(num2str(i1))';

end

function swc = buildSWC(tree)

% builds swc from tree format
nN = length(tree.X);
swc = zeros(nN,7);
swc(:,1) = 1:nN;
swc(:,3) = tree.X;
swc(:,4) = tree.Y;
swc(:,5) = tree.Z;
swc(:,6) = tree.D / 2; % radius

PIdx = find(tree.dA);
[parent, child] = ind2sub(size(tree.dA),PIdx);
P = zeros(nN,1);
P(parent) = child;
P(1) = -1;
swc(:,7) = P;
swc(1,2) = 1;

% set type based on number of connections
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


function filelist = loadfilelist(filename)
% filelist = loadfilelist(filename)
% Read a plain text file for all image names. One line is an image name.
%
% By Hanchuan Peng
% Jan,2001
% June, 2005. Fix the non-return value bug

filelist = [];
fid = fopen(filename);
if fid==-1,
    disp(['Error to open the file : ' filename]);
    return;
else
    i=1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break; end;
        filelist{i} = deblank(tline);
        i = i+1;
    end;
end;
fclose(fid);
end


function a = load_v3d_neuron_file(filename, b_minusFirst)
%function a = load_v3d_neuron_file(filename, b_minusFirst)
% Load the swc file as a neuron structure
% by Hanchuan Peng
% 070726
% 070809: I notice the majority of neuron swc files have their soma have
% the corredinates (0,0,0), but there are still exceptions. Thus I subtract
% the soma coordinate to assure any neurons have the same origin.
% 080513: add a b_minusFirst flag
% 0904: add a get first 7 columns check
% 090415: add a deblank check

if nargin<2,
    b_minusFirst=0;
end;

L = loadfilelist(filename);
a = zeros(length(L), 7);

k=0;
for i=1:length(L),
    if isempty(deblank(L{i})),
        continue;
    end;
    if (L{i}(1)=='#'),
        continue;
    end;
    
    k=k+1;
    tmp = str2num(L{i});
    a(k,:) = tmp(1:7);
end;

a = a(1:k,:); % remove the non-used lines

% make sure all the origin (neuron soma) will be 0
if b_minusFirst,
    a(:,3:5) = a(:,3:5) - repmat(a(1,3:5), size(a,1), 1);
end;

return;
end