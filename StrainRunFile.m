% This file takes a swc file, smooths it, and returns the 
% local/global strains for each neuron.

% If used please cite:
% E. Bar-Kochba, M. Scimone, J. Estrada, and C. Franck, "Strain and 
% rate-dependent diffuse axonal injury of 3d neuron cultures under 
% compression," Biophysical Journal, vol. 3, no. 110, p. 320a, 2016.

clear, clc, close all

% Add all m-files to search path
addpath(genpath(cd))
F2 = uigetdir(cd,'Select where swc files are stored files');
cd(F2);
swcfiledir = dir('*.swc');
swcfilename = swcfiledir.name;
% Farfield strain
nu = 0.056; % Poisson's Ratio
mag = 0.3;  % Strain magnitude 
Einf = mag*[-1*nu      0         0; 
                0      -1*nu     0; 
                0      0         1];
% Change signs depending on frame of reference and known input.
% Strains


[E_local, E_mean, K_mean, Ef, tree, Ef_eig] = strainsTree(swcfilename,Einf);

% Before plotting with either option, one can smooth node diameters.
% (Not required) Editing the model in neuron studio is another option. 
for j = 1:length(tree)
m = prctile(tree{j}.D,90);
tree{j}.D(tree{j}.D < m) = m;
rf = findSeg({tree{j}}); rf = rf{1};
close all
for ii = 1:size(rf,1) 
    seg = nonzeros(rf(ii,:));
    tree{j}.D(seg) = smooth(tree{j}.D(seg));
end
end

% Plot strains

neuron = 1; % first neuron in tree structure
% component: which strain component to plot
%               1 - Radial Strain
%               2 - Radial Shear Strain
%               3 - Axial Strain
%               4 - Axial Shear Strain
component = 3; % axial strain
plotStrain(E_local{j},tree{j},neuron,component)


% Plot strains using Tecplot Focus

tic
neuron = 1; % first neuron in tree structure
component = 3; % axial strain
savename = 'networkreduce.plt'; % plt file name
plotTec(E_local{j},tree{j},neuron,component,savename)
toc


%% Finding strain values in region of interest

% This block is meant to be run with the network example data included.
% Preceding block must be run first. 

% Use segmented tiff
I0 = loadtiff('140714_impact_xy03_3_Segmented3.tif');

% Micrometer to voxel conversion factor
um2vxl = [0.0776851 0.0776851 0.5];

% Show selected ROI
plot = 1;

% ROI strains
[ROI_local,ROI_global] = strainROI(I0,um2vxl,tree,E_local,Ef_eig,plot);
