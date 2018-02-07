function plotTec(E_local,tree,neuron,component,savename)
%
% plotTec(E_local,tree,neuron,component)
%
% Create a plt file ready for use in Tecplot focus. The strain component of
% the selected model neuron will be plotted. 
%
% INPUTS
%--------
% E_local:   local strain values 
% tree:      tree structure defining model neuron
% neuron:    which neuron to plot if tree is a cell array
% component: which strain component to plot
%               1 - Radial Strain
%               2 - Radial Shear Strain
%               3 - Axial Strain
%               4 - Axial Shear Strain
% savename:  plt file name

% Take arrays out of cells
if iscell(tree), tree = tree{neuron}; end
if iscell(E_local), E_local = E_local{neuron}; end

% Set up tdata structure. See mat2tecplot for detailed documentation. 
tdata = [];
tdata.Nvars = 4;
tdata.vformat(1:tdata.Nvars) = 2;
switch component
    case 1, tdata.varnames = {'x','y','z','Radial Strain'};
    case 2, tdata.varnames = {'x','y','z','Radial Shear Strain'};
    case 3, tdata.varnames = {'x','y','z','Axial Strain'};
    case 4, tdata.varnames = {'x','y','z','Axial Shear Strain'};
end

% Set soma NaN value to max strain value. Branch and endpoints NaN values
% are not shown in plot. See notice in strainsTree.m for why NaN values
% exist.
E_local(1,:) = max(E_local(:,component));

% All other NaN values are set to the same strain value as their parent
% node
% Find parent nodes (see buildSWC) 
PIdx = find(tree.dA);
[parent, child] = ind2sub(size(tree.dA),PIdx);
P = zeros(length(tree.X),1);
P(parent) = child;
P(1) = -1; % set soma to -1
% Substitute parent node values
k = 0;
while sum(isnan(E_local(:,component))) > 0
    idnan = find(isnan(E_local(:,component)));
    E_local(idnan,component) = E_local(P(idnan),component);
    % Saftey in case while loop takes too long 
    k = k + 1; % Counter
    if k == 100
        E_local(idnan,component) = median(E_local(:,component));
        break
    end
end

% Create a sphere
[xu,yu,zu] = sphere;
for ii = 1:length(tree.X)
    
    % Modify sphere for each node
    x = xu*tree.D(ii)/2 + tree.X(ii);
    y = yu*tree.D(ii)/2 + tree.Y(ii);
    z = zu*tree.D(ii)/2 + tree.Z(ii);
    
    [~,nv] = reducepatch(surf2patch(x,y,z),0.2);
    x = nv(:,1); y = nv(:,2); z = nv(:,3);
    tri = delaunay(nv);
    
    % Save volume information in tdata structure
    tdata.FEvolumes(ii).x = x';
    tdata.FEvolumes(ii).y = y';
    tdata.FEvolumes(ii).z = z';
    tdata.FEvolumes(ii).v(1,:) = ones(size(z))*E_local(ii,component);
    tdata.FEvolumes(ii).e2n = tri;
    tdata.FEvolumes(ii).zonename = num2str(ii);
        
end

mat2tecplot(tdata,savename)