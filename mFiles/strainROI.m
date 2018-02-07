function [ROI_local,ROI_global] = strainROI(I0,um2vxl,tree,E_local,Ef_eig,plot)
%
% strainROI(I0,um2vxl,E_local,Ef_eig)
%
% Determine the local and global strains at nodes within the user-defined
% ROI. Is meant for a network of neurons, but can be used for single model
% neurons or anytime tree input is a structure, not a cell array.
%
% INPUTS
%------------------
% I0:           Segmented volumetric image
% um2vxl:       Micrometer to voxel conversion factors in x,y,z
% tree:         Tree structure as a structure class
% E_local:      Local strains as an array
% Ef_eig:       Local strains as an array but calculated at the edges
% plot:         Logical input to determine whether or not to plot
%--------          
% OUTPUTS            
% -----------------
% ROI_local:    Local strains of nodes within ROI as an array in the same
%               format as E_local
% ROI_global:   Mean strain of entire ROI
%

imshow(max(I0,[],3));
% Need to switch x y because of how ginputc works
[y,x] = ginputc(2,'Color','r');
close all
x = sort(x);
y = sort(y);

% Show selected ROI
if plot
RGB = insertShape(max(I0,[],3),'FilledRectangle',[y(1) x(1)...
    y(2)-y(1) x(2)-x(1)]);
imshow(RGB)
end

% Define region of interest
x = x*um2vxl(1);
y = y*um2vxl(2);

% Find nodes within region of interest
idx = find(tree.X > x(1) & tree.X < x(2));
idy = find(tree.Y > y(1) & tree.X < y(2));
C = intersect(idx,idy); % Indices that meet both x and y region of interest

% Pull out the local strains
ROI_local = E_local(C,:);

% Segment lengths
len = len_tree(tree);
len = len(C);                      
L = nansum(len);

% Mean strain
Ef_eig = Ef_eig(C,:);
ROI_global = zeros(1,size(Ef_eig,2));
for i = 1:size(Ef_eig,2)
    ROI_global(1,i) = nansum(Ef_eig(:,i).*len)/L;
end