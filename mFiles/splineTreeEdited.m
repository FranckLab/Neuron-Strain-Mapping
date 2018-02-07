function [t,n,b,K,T] = splineTreeEdited(swc,location,plot01)
%
% splineTree(swc,method,location,plot01)
%
% Performs cubic spline interpolation on the each segment of the neuron, 
% and calculates Frenet Frame at each node
%
% INPUTS
%------------------
% swc:          Matrix in swc format rather than tree structure format.
%               Output from local function buildSWC is in the proper format
%
% location:     Vertex 'v' (node) or edge 'e' (along the segment)
%
% plot01:       Logical value whether or not to plot
%--------              
% OUTPUTS             
% -----------------   
% t:            Tangent unit vector at each node. n x 3 matrix where n is
%                   the number of nodes
% n:            Normal unit vector at each node
% b:            Binormal unit vector at each node
% K:            Curvature at each node
% T:            Torsion at each node
% 

if nargin < 3, plot01 = 0; end
if nargin < 2, location = 'v'; end
if isempty(location), location = 'v'; end

% Only frenet frame is implemented
method = 'frenet';

[t,n,b,K,T] = frenetTree_fun(swc,method, location,plot01);

end



function [t,n,b,K,T] = frenetTree_fun(swc,~, location,plot01)

% Get all branches
BR = branchTree(swc); 
nBR = length(BR);

% Cartesian coordinates of nodes
r = swc(:,3:5);

% Make NaN
t = zeros(size(r,1),3);
n = t;
b = t;
K = t(:,1);
T = t(:,1);

for i = 1:nBR % number of branches
    
    % Get current branch
    BRi = BR{i};
    
    % Get coordinates of nodes in current branch
    ri = r(BR{i},:);
    
    % Don't pass branches less than three elements or else error in
    % derivatives later
    if size(ri,1) > 3
        [tangent, normal, binormal, kappa, tau, xpp, pp] = ...
            splineFrenet(ri,location);

        if plot01
            figure(99);
            xpp_ = linspace(xpp(1),xpp(end),10000);
            ripp = ppval(pp,xpp_)';
            hold on;
            plot3(ri(:,1),ri(:,2),ri(:,3),'c.');
            plot3(ripp(:,1),ripp(:,2),ripp(:,3),'k-','linewidth',4);
            axis image
            view(3)
            grid on;
            box on;           
            quiver3(ri(:,1),ri(:,2),ri(:,3),tangent(:,1),tangent(:,2),...
                tangent(:,3),'r')
            quiver3(ri(:,1),ri(:,2),ri(:,3),normal(:,1),normal(:,2),...
                normal(:,3),'g')
            quiver3(ri(:,1),ri(:,2),ri(:,3),binormal(:,1),binormal(:,2),...
                binormal(:,3),'b')
            legend('Nodes','Spline','Tangent','Normal','Binormal')
        end

        if strcmpi(location(1),'e')
            BRi = BRi(2:end);
            t(BRi,:) = tangent;
            n(BRi,:) = normal;
            b(BRi,:) = binormal;
            K(BRi) = kappa;
            T(BRi) = tau;
        else
            % valid vertex locations 
            t(BRi,:) = tangent;
            n(BRi,:) = normal;
            b(BRi,:) = binormal;
            K(BRi) = kappa;
            T(BRi) = tau;
        end
    end
end
end




function branch = branchTree(swc)

% Integer representing the type of neuronal segment, such as soma, axon...
T = swc(:,2);

% Parent (the integer label) of the current point or -1 to indicate an 
% origin (soma).
P = swc(:,7);

% start at end point (6) or fork (5)
start = find(T == 6 | T == 5 | T == 0);
nStart = length(start);

branch = cell(nStart,1);
for i = 1:nStart
    
    branch{i}(1,1) = start(i);
    branch{i}(2,1) = P(start(i));
    while T(branch{i}(end,1)) == 3 % continue while node is dendrite (T = 3)
        branch{i}(end + 1,1) = P(branch{i}(end));
    end
    branch{i} = flipud(branch{i});  % direct graph away from root
end
end
