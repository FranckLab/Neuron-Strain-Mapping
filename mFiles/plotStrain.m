function plotStrain(E_local,tree,neuron,component)
%
% plotStrain(E_local,tree,neuron,component)
%
% Plot a strain component of a model neuron
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
% 

% Take arrays out of cells
if iscell(tree), tree = tree{neuron}; end
if iscell(E_local), E_local = E_local{neuron}; end

% Set up figure window
figure
set(gcf, 'Color', 'w');
hold on;

% Set soma NaN value to max strain value. Branch and endpoints NaN values
% are not shown in plot. See notice in strains tree for why NaN values
% exist.
E_local(1,:) = max(E_local(:,component));

% Plot spheres with strain values as cdata in surf()
[xu,yu,zu] = sphere;
for ii = 1:length(tree.X)
  x = xu*tree.D(ii)/2 + tree.X(ii);
  y = yu*tree.D(ii)/2 + tree.Y(ii);
  z = zu*tree.D(ii)/2 + tree.Z(ii);
  c = ones(size(z))*E_local(ii,component);
  surf(x,y,z,c,'edgecolor','none');
end

% Make the correct title
switch component
    case 1, title('Radial Strain','fontsize',16,'fontweight','bold')
    case 2, title('Radial Shear Strain','fontsize',16,'fontweight','bold')
    case 3, title('Axial Strain','fontsize',16,'fontweight','bold')
    case 4, title('Axial Shear Strain','fontsize',16,'fontweight','bold')
end

% Figure options
h = colorbar;
h.LineWidth = 2;
h.FontSize = 12;
h.FontWeight = 'bold';
view(3);
box on
ax = gca;
ax.BoxStyle = 'full';
set(0,'DefaultFigureColor','remove')
set(gca,'xtick',[],'ytick',[],'ztick',[])
axis equal;