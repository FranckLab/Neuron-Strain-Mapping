function volPP = PreProcessNeuron3D(vol, KernelSize, Iterations,...
    PreFix, SaveMat, SaveTif, Plot)
% Syntax: volPP = PreProcessNeuron3D(vol, SavePrefix, KernelSize, SaveMat, SaveTif, Plot)
% INPUTS:
%  1) vol: volume
%  2) SaveMat: save the variable as a .mat file? (Yes = 1)
%  3) SaveTif: save the variable as a .tif file? (Yes = 1)
%  4) Plot: plot a maximum Z-projection of the preprocessed volume
%
%  Author: Eyal Bar-Kochba, Brown University (2012) <><

if nargin < 6, Plot = 0; end
if nargin < 5, SaveTif = 0; end
if nargin < 4, SaveMat = 0; end
if nargin < 3, Iterations = 1; end


volPP = vol - mean2(vol); %% Subtract Background
for i = 1:Iterations
    if Iterations == 0, break; end
    volPP = medfilt3(volPP,KernelSize); %% Apply median filter
end
%volPP = imadjust3(volPP); %% Adjust Contrast THIS IS BINARIZING THE IMAGE

if Plot == 1
% Plot maximum Z-projection for comparison
figure
subplot(1,2,1); imagesc(squeeze((max(vol,[],3)))), colormap gray, axis image
subplot(1,2,2); imagesc(squeeze((max(volPP,[],3)))), colormap gray, axis image
end

if SaveMat == 1, save(PreFix ,'volPP'), end
if SaveTif == 1, mat2tiff(volPP,'xy',class(vol),PreFix), end
end

% Contrast adjustment
function vol = imadjust3(vol)

ContrastLimits = stretchlim(max(vol,[],3));
for k = 1:size(vol,3), vol(:,:,k) = imadjust(vol(:,:,k),ContrastLimits,[]); end
end

% Median filter
% function vol = medfilt3(vol, KernelSize)
% KernelSize = ones(1,2)*KernelSize;
% for k = 1:size(vol,3), vol(:,:,k) = medfilt2(vol(:,:,k),KernelSize); end
% end

    