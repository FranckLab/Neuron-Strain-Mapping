% This file takes nd2 files and returns processed tiff files ready to
% import into Neuron Studio. This file is meant to be run block by block
% as needed by the user. 

% If used please cite:
% E. Bar-Kochba, M. Scimone, J. Estrada, and C. Franck, "Strain and 
% rate-dependent diffuse axonal injury of 3d neuron cultures under 
% compression," Biophysical Journal, vol. 3, no. 110, p. 320a, 2016.

clear, clc, close all

% Add all m-files to search path - edit accordingly
addpath(genpath('E:\Documents\Franck lab\TBIhypotherm\NeuronSegmentationTestPackage'))

%% Convert nd2 files to mat 

% Prefix of nd2 files to be converted
%filename = '140729_Impact*';

% Select folder to save mat files to
%Make sure there are no subdirectories in this folder!
F1 = uigetdir(cd,'Select folder to save mat files to');
savePath = F1;

% True/False statement to use parfor loop
nProcessors = 1;

% Select folder where nd2 files are stored
oldfolder = cd;
F1 = uigetdir(cd,'Select folder where nd2 files are stored');
cd(F1)

% Convert to mat files
fileDir = dir('*.nd2');
filename = {fileDir.name};
%Sometimes bfopen organizes zsteps such that step alternates by channel 
%and not z-position. For this reason we recommend separating channels beforehand.   
    for j = 1:length(filename)
        data{j} = bfopen(filename{j});
    end
    
xpx = length(data{1, 1}{1, 1}{1, 1});
ypx = length(data{1, 1}{1, 1}{1, 1});
zslices = length(data{1, 1}{1, 1});
I = zeros (xpx,ypx,zslices);


for ii = 1:length(data)
    DataFolder = filename{ii}(1:end-4);
    mkdir(savePath,DataFolder)
    newSavePath = num2str([savePath,filesep,DataFolder]);
    for jj = 1:size(data{1},1)
        dir_to_save = [newSavePath,filesep,filename{ii}(1:end-4),'_',num2str(jj)];
        count=1;
        for kk = 1:size(data{1,1}{1,1},1)
            I(:,:,count) = data{ii}{jj, 1}{kk, 1};
            
            count = count+1;
        end
        save(dir_to_save,'I');
    end
    cd(savePath)
end

%bioformats2mat([filename,'.nd2'], savePath, nProcessors);

% cd to original dir
cd(oldfolder)

%% Load mat files for pre-processing

% This block loads one folder of mat files at a time into variable mydata

% Select where mat files are saved
oldfolder = cd;
F1 = uigetdir(cd,'Select folder where mat files are saved');
cd(F1)

% This uses fileparts and dir to find out how many files are in folder

% NOTE: Assumes filename prefix is the same as folder name (if using
% bioformats2mat, this should automatically be done)

% Example: folder = AnnulusNice40x 
%          file   = AnnulusNice40x_t01.mat

% If they're not the same, then comment out fileparts and type prefix of
% files as a string into myfolder in dir function

[~,myfolder,~] = fileparts(cd);
myfolderinfo = dir([myfolder,'*.mat']);
numfiles = length(myfolderinfo);
filenames = {myfolderinfo.name}.';

% This stores volume stacks of each time point in a cell
% of size numfiles x 2 into variable mydata
mydata = cell(numfiles,1);
for ii = 1:numfiles
    namefile = myfolderinfo(ii).name;
    vol = load(namefile);
    names = fieldnames(vol);
    fieldname = names{1};
    if length(vol.(fieldname)) == 2        
        for jj = 1:2
            mydata{ii,jj} = vol.vol{jj};      
        end       
    else
        mydata{ii} = vol.I;   
    end
end

cd(oldfolder)
%still full bitdepth
clearvars -except mydata filenames


%% Correct for drift

% This block is optional for use cases where the researcher suspects there
% is drift in the system over time. Volume stacks specified by filename
% is used to correct for experiment drift by using the median values of
% u1 u2 u3 outputs from Fast Iterative Digital Volume Correlation(FIDVC).
% Variable mydata is translated accordingly. 

% Select where mat file is saved
addpath(genpath(cd))
oldfolder = cd;
F1 = uigetdir(cd,'Select folder where mat files are saved');
cd(F1)

% FIDVC inputs
sSize = [128 128 64];
incORcum = 'incremental';
filename{1} = '140729_Impact01_t01.mat';
filename{2} = 1;

% Estimate displacements via IDVC
[u,~,~] = funIDVC(filename, sSize, incORcum);

cd(oldfolder)

% Translate volume stack according to median values of u1 u2 u3
D = zeros([size(mydata{1}),3]);
for ii = 1:length(u)
   
    % Find median values of u1 u2 u3
    D(:,:,:,1) = median(u{ii}{1,1}(:));
    D(:,:,:,2) = median(u{ii}{1,2}(:));
    D(:,:,:,3) = median(u{ii}{1,3}(:));  
    
    % Translate both live and dead channel
    mydata{ii+1,1} = imwarp(mydata{ii+1,1},D);
    mydata{ii+1,2} = imwarp(mydata{ii+1,2},D);
    
end

%% Image processing

% This block filters noise using a median filter, applies gamma correction, 
% and thresholds the background of each volume stack. Filtered data is
% stored in variable mydatafilt.

% Which channel contains signal of interest 
channel = 1;
% Number of times median filter is applied
Iterations = 2;
% Neighborhood size of median filter
KernelSize = [3 3 3];
mydatafilt = cell(size(mydata,1),1);
for ii = 1:size(mydata,1)
    mydatafilt{ii,1} = PreProcessNeuron3D(mydata{ii,channel},...
        KernelSize,Iterations);
end

% View maxiumum intensity projection of the first volume
% imshow(max(mydatafilt{1},[],3))

%% Convert from mat to tiff files

% This block converts variable mydatafilt into tiff files one channel at a 
% time.

% Select where to save tiff files. You may want to make a new folder to
% keep files organized
oldfolder = cd;
F1 = uigetdir(cd,'Select where to save tiff files');
cd(F1)

for ii = 1:length(mydata)

    % Specify volume stack and coonvert bitdepth
    vol = uint8(mydatafilt{ii,1});

    % Direction to save volume. See mat2tiff.m for more options. 
    saveslice = 'xy';

    % Specify numeric class
    numclass = 'uint8';

    % Specify filename
    filename = filenames{ii}(1:end-4);

    mat2tiff(vol,saveslice,numclass, filename)
end

cd(oldfolder)