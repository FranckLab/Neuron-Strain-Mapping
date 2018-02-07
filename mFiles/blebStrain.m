function [bleb_local,bleb_global,C] = ...
    blebStrain(treeBleb,tree,E_local,Ef_eig,plot01,component)
%
% blebStrain(treeBleb,tree,E_local,Ef_eig,plot)
%
% Determine the local and global strains of tagged blebs
%
% INPUTS
%------------------
% treeBleb:         Tree structure of individual blebs. Can be a class
%                       structure or cell.
% tree:             Tree structure of entire model neuron. Must be a class
%                       stucture.
% E_local:          Local strains for entire model neuron
% Ef_eig:           Local strains for entire model neuron calculated along
%                       segment edges
% plot:             Logical 1 x 2 vector determining whether to plot mean
%                       bleb strain values and local bleb strain values,
%                       respectively
% component:        Which strain component to plot
%                        1 - Radial Strain
%                        2 - Radial Shear Strain
%                        3 - Axial Strain
%                        4 - Axial Shear Strain
%--------          
% OUTPUTS            
% -----------------
% bleb_local:       Strain values at nodes that make up each bleb. The
%                       nodes of each bleb are stored in a cell. Within
%                       each cell, the strains are in the same format as
%                       E_local.
% bleb_global:      Mean strain values at each bleb. Stored as an m x n
%                       array where m is the number of blebs and the 
%                       columns represent strain components
% C:                Indices of nodes that were used in bleb strain
%                       calculations. Stored as a cell array. 

% Search radius
for ii = 1:size(treeBleb,1), r(ii) = treeBleb{ii}.D; end
radius = mean(r)/2;

% Pre-Allocate
bleb_local = cell(size(treeBleb,1),1);

% Find bleb strains
len = len_tree(tree);

for ii = 1:size(treeBleb,1)
    
    % Range of values to search
    x = [treeBleb{ii}.X - radius treeBleb{ii}.X + radius];
    y = [treeBleb{ii}.Y - radius treeBleb{ii}.Y + radius]; 
    z = [treeBleb{ii}.Z - radius treeBleb{ii}.Z + radius];
    
    % Find which nodes in x,y, OR z fit match values
    idx = find(tree.X > x(1) & tree.X < x(2));
    idy = find(tree.Y > y(1) & tree.Y < y(2));
    idz = find(tree.Z > z(1) & tree.Z < z(2));
    
    % Find matching coordinates across x,y, AND z
    C{ii} = intersect(intersect(idx,idy),idz);
    
    % Grab the strain values at these nodes
    bleb_local{ii} = E_local(C{ii},:);
    
    % Mean strain at each bleb point
    lenbleb = len(C{ii});
    L = nansum(lenbleb);
    Ef_eigbleb = Ef_eig(C{ii},:);
    for i = 1:size(Ef_eig,2)
        bleb_global(ii,i) = nansum(Ef_eigbleb(:,i).*lenbleb)/L;
    end
  
end

% Plot mean strain at blebs
if plot01(1)
figure,
set(gcf,'color','w');
scatter(1:size(bleb_global,1),bleb_global(:,component),75,'filled')
switch component
    case 1, title('Bleb Radial Strain','fontsize',16,'fontweight','bold')
    case 2, title('Bleb Radial Shear Strain','fontsize',16,'fontweight','bold')
    case 3, title('Bleb Axial Strain','fontsize',16,'fontweight','bold')
    case 4, title('Bleb Axial Shear Strain','fontsize',16,'fontweight','bold')
end
xlabel('Bleb Index Number','FontSize',14,'FontWeight','bold')
ylabel('Strain','FontSize',14,'FontWeight','bold')
set(gca,'xtick',1:length(treeBleb))
xlim([0 length(treeBleb)+1])
set(gca,'fontWeight','bold','fontsize',12)
box on 
end

% Plot strain at nodes that make up blebs
if plot01(2)
figure,
set(gcf,'color','w');
hold on
for ii = 1:size(bleb_local,1) 
    scatter(1:length(bleb_local{ii}(:,component)),...
        bleb_local{ii}(:,component),75,'filled')
end
switch component
    case 1, title('Nodal Radial Strain','fontsize',16,'fontweight','bold')
    case 2, title('Nodal Radial Shear Strain','fontsize',16,'fontweight','bold')
    case 3, title('Nodal Axial Strain','fontsize',16,'fontweight','bold')
    case 4, title('Nodal Axial Shear Strain','fontsize',16,'fontweight','bold')
end
xlabel('Nodal Index Number','FontSize',14,'FontWeight','bold')
ylabel('Strain','FontSize',14,'FontWeight','bold')
set(gca,'fontWeight','bold','fontsize',12)
[nrows, ~] = cellfun(@size, bleb_local);
set(gca,'xtick',1:max(nrows))
box on
% [h, ~] = legend('Bleb 1','Bleb 2','Bleb 3');
end

end