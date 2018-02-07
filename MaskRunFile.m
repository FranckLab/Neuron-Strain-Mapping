% This file creates a mask from a swc file and filters out all none
% model nuerons.

% If used please cite:
% E. Bar-Kochba, M. Scimone, J. Estrada, and C. Franck, "Strain and 
% rate-dependent diffuse axonal injury of 3d neuron cultures under 
% compression," Biophysical Journal, vol. 3, no. 110, p. 320a, 2016.

clear, clc, close all

% Adjust search path
addpath(genpath(cd))

% Select folder where tiffs to be segmented are stored. Adjust string name
% as needed.
F1 = uigetdir(cd,'Select where processed tiffs are stored');
cd(F1);
filedir = dir('*.tif');

% Specify swc file path
F2 = uigetdir(cd,'Select where swc files are stored files');
cd(F2);
swcfiledir = dir('*.swc');
swcfilename = swcfiledir.name;


% Specify where to save segmented tiffs
F3 = uigetdir(cd,'Specify where to save segmented tiffs');
cd(F3);
SegTifdir = dir('*.swc');

cd(F1)

% Specify micrometer to voxel conversion factor
um2vxl = [0.173 0.173 1];  % 0.0776851 0.1553701

% Specify how many voxels the mask will be dilated by
vDilate = 10;

for ii = 1:length(filedir) % number of tiff files
    
% Load processed tiff file
I0 = loadtiff(filedir(ii).name);

if ii == 1
    % Create a mask
    BW = swc2mask(swcfilename, size(I0), um2vxl);
    % Find indices in BW within distance d 
    D = bwdist(BW,'euclidean');
    ind = find(D <= vDilate);
    % Dilate mask accordingly
    BW(ind) = 1;
    % Check numeric class of I0
    if isa(I0,'uint8') == 1
        BW = im2uint8(BW)/255;
    else
        BW = uint16(BW);
    end
end

% Apply mask
Seg = I0.*BW;

% Save segmented tiff file
cd(F3)
savename = [swcfilename(1:end-4),'_',num2str(ii)];
mat2tiff(Seg,'xy','uint8',savename);
cd(F1)
end

% Maximum Intesnity Project along 3 axes
% MIPz = max(Seg, [], 3);
% MIPy = squeeze(max(Seg, [], 2));
% MIPx = squeeze(max(Seg, [], 1));
% convertion = um2vxl(3)/um2vxl(1);
% [newmeshx, newmeshy]= meshgrid (1:1024, 1:1/convertion:77);
% newMIPy = interp2(double(MIPy)', newmeshx,newmeshy);
% newMIPx = interp2(double(MIPx)', newmeshx,newmeshy);
% 
% imshow(MIPy,[])