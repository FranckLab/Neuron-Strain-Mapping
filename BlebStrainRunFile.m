% This file finds the axial strain, axial shear strain, radial
% strain, radial shear strain, and curvature of all nodes within a distance 
% of the mean bleb radius.

% If used please cite:
% E. Bar-Kochba, M. Scimone, J. Estrada, and C. Franck, "Strain and 
% rate-dependent diffuse axonal injury of 3d neuron cultures under 
% compression," Biophysical Journal, vol. 3, no. 110, p. 320a, 2016.

clear, clc, close all

% Add all m-files to search path
addpath(genpath(cd))

% Farfield strain
nu = 0.056; % Poisson's Ratio
mag = 0.3;  % Strain magnitude
Einf = mag*[-1*nu      0         0; 
                0      -1*nu     0; 
                0      0         1];

% Strains across entire neuron
[E_local, ~, ~, ~, tree, Ef_eig] = strainsTree('40x.swc',Einf);

% Read in swc file that marks blebs
[~,treeBleb] = readSWC('tag3blebs.swc');

% Plot mean (global) axial bleb and local bleb strain values
plot01 = [1 1];

% Axial strain component
component = 3;

% Bleb strains
[bleb_local,bleb_global,C] = ...
    blebStrain(treeBleb,tree,E_local,Ef_eig,plot01,component);

% Show which nodes were selected during calculations
%{
um2vxl = [0.1553701 0.1553701 0.5];
I0 = loadtiff('40xSegmented.tif');
imshow(max(I0,[],3))
hold on
for ii = 1:length(C)
plot(tree.X(C{ii})/um2vxl(1),tree.Y(C{ii})/um2vxl(2),'r.','MarkerSize',20)
end
%}