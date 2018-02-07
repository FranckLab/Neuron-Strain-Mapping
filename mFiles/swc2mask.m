function BW = swc2mask(swc, sizeI, um2vxl)
%
% swc2mask(swc,sizeI,um2vxl)
%
% Create a logical mask of sizeI from the swc file and scale it according to
% um2vxl
%
% INPUTS
%------------------
% swc:          swc file from Neuron Studio
% sizeI:        size of volumetric image the swc file is associate with
% um2vxl:       1x3 vector specifying the x,y,z conversion factors
%--------
% OUTPUTS
% -----------------
% BW:           logical mask of sizeI created from swc file 
%

% Read in tree structure
swc = readSWC(swc);

% Make sure swc is a cell
if ~iscell(swc), swc = {swc}; end

% Number of model neurons in swc file
nSoma = length(swc);

% Create mask
BW = false(sizeI);
for i = 1:nSoma
    % Create mask for one neuron at a time
    BW_ = swc2mask_call(swc{i},BW,um2vxl);
    % Keep all model neurons using | operator
    BW = BW_ | BW;
end

end

% Create logical mask of one neuron 
function BW = swc2mask_call(swc, BW, um2vxl)

% Extract swc parts. Read readSWC.m documentation for more details on what
% each column represents.
N = swc(:,1);
R = swc(:,6);
P = swc(:,7);
X = ceil([swc(:,4)/um2vxl(2) swc(:,3)/um2vxl(1) swc(:,5)/um2vxl(3)]);
R = ceil(R/um2vxl(1));

% Start with soma. Creates a logical vector length(P), with true at value
% = -1 (soma index)
somaIdx = P == -1;

% Indices of upper left corner of sphere representing soma
OS = X(somaIdx,:)-R(somaIdx);

% Create a logical sphere of soma
hull = fillSphere(R(somaIdx), um2vxl);

% Indices that represent diameters of the sphere in every dimension
xIdx = (1:size(hull,2)) + OS(2);
yIdx = (1:size(hull,1)) + OS(1);
zIdx = (1:size(hull,3)) + OS(3);

% Size of tiff file
sizeI = size(BW);

% Remove indices of the sphere that are not captured by the original image
badIdx = xIdx < 1 | xIdx > sizeI(2);
xIdx(badIdx) = [];
hull(:,badIdx,:) = [];

% Repeat for y
badIdx = yIdx < 1 | yIdx > sizeI(1);
yIdx(badIdx) = [];
hull(badIdx,:,:) = [];

% Repeat for z
badIdx = zIdx < 1 | zIdx > sizeI(3);
zIdx(badIdx) = [];
hull(:,:,badIdx) = [];

% Place the soma into binary BW array
BW(yIdx,xIdx,zIdx) = hull;

fprintf(1,['swc2mask (%):']);
for n = 2:length(N) % number of nodes in current neuron
    
    percentage = ceil(n/length(N)*100);
    fprintf(1,['\b\b\b\b\b%3.2i%%\n'],percentage);
    
    % Coordinates and radius of the next sphere and its parent 
    X1 = X(n,:); R1 = R(n);
    X2 = X(P(n),:); R2 = R(P(n));
    
    % Create convex hull from the two spheres
    hull = ConvexHullSphere(X1,R1,X2,R2,um2vxl);
    
    % Remove of the convex hull not captured by the original image
    OS = min(X1-R1,X2-R2) - 1;
    xIdx = (1:size(hull,2)) + OS(2);
    yIdx = (1:size(hull,1)) + OS(1);
    zIdx = (1:size(hull,3)) + OS(3);
    
    badIdx = xIdx < 1 | xIdx > sizeI(2);
    xIdx(badIdx) = [];
    hull(:,badIdx,:) = [];
    
    badIdx = yIdx < 1 | yIdx > sizeI(1);
    yIdx(badIdx) = [];
    hull(badIdx,:,:) = [];
    
    badIdx = zIdx < 1 | zIdx > sizeI(3);
    zIdx(badIdx) = [];
    hull(:,:,badIdx) = [];
    
    % Place conv hull into BW
    BW(yIdx,xIdx,zIdx) = BW(yIdx,xIdx,zIdx) | hull;
    
end
end

% Create a logical sphere
function BW = fillSphere(R,um2vxl)

scale = um2vxl(3)/um2vxl(1);
idx = -R:R;
[mx, my, mz] = meshgrid(idx,idx,idx);
mz = mz*scale;
BW = sqrt((mx).^2 + (my).^2 + (mz).^2) <= R;

end