function [tangent, normal, binormal, kappa, tau, t, s] = splineFrenet(x,location)
%
% splineFrenet(x,location)
%
% Perform cubic spline interpolation in a set of 3D points and find the
% corresponding Frenet frame at each of those points. 
%
% INPUTS
%------------------
% x:            n x 3 break points for spline where n is the number of
%               points and the three columns are the cartesian coordinates
%
% location:     Vertex 'v' (node) or edge 'e' (along the segment)
%--------          
% OUTPUTS            
% -----------------
% tangent:      Tangent unit vector at each node. n x 3 matrix where n is
%                   the number of nodes
% normal:       Normal unit vector at each node
% binormal:     Binormal unit vector at each node
% kappa:        Curvature at each node
% tau:          Torsion at each node
% t:            query points along arclength
% s:            spline function
%

if nargin < 2, location = 'v'; end

ts = [0; cumsum(sqrt(sum(diff(x,[],1).^2,2)))]; % arclength

% Cubic spline interpolation
s = spline(ts',reshape(x',[3, length(x)]));

% Modify arclength if edges are chosen
if strcmpi(location(1),'e')
   t = diff(ts)/2 + ts(1:end-1);  % midpoints of edges
else
   t = ts;
end

if size(s.coefs,1) == 3
    s.coefs = s.coefs';
end

% first derivative, pp form
s1 = s;
s1.order = s.order-1; 
s1.coefs = bsxfun(@times, s.coefs(:,1:end-1), s1.order:-1:1);

% second derivative, pp form
s2 = s1;
s2.order = s1.order-1; 
s2.coefs = bsxfun(@times, s1.coefs(:,1:end-1), s2.order:-1:1);

% third derivative, pp form
s3 = s2;
s3.order = s2.order-1; 
s3.coefs = bsxfun(@times, s2.coefs(:,1:end-1), s3.order:-1:1);

dr = ppval(s1,t); % first derivative at x
drr = ppval(s2,t); % second derivative at x
drrr = ppval(s3,t); % third derivative at x

dx = dr(1,:)';
dy = dr(2,:)';
dz = dr(3,:)';

dxx = drr(1,:)';
dyy = drr(2,:)';
dzz = drr(3,:)';

dxxx = drrr(1,:)';
dyyy = drrr(2,:)';
dzzz = drrr(3,:)';

% r' x r''
drxdrr1 = (dzz.*dy - dyy.*dz);
drxdrr2 = (dxx.*dz - dzz.*dx);
drxdrr3 = (dyy.*dx - dxx.*dy);

% calculate curvature kappa = |r' x r''| / |r'|^3
num = sqrt( drxdrr1.^2 + drxdrr2.^2 + drxdrr3.^2 );
den = (dx.^2 + dy.^2 + dz.^2).^(3/2);
kappa = num./den;
kappa = kappa(:);

% calculate torsion tau = ((r' x r) dot r''') / |r' x r''|^2
num = dxxx.*drxdrr1 + dyyy.*drxdrr2 + dzzz.*drxdrr3;
den = drxdrr1.^2 + drxdrr2.^2 + drxdrr3.^2;
tau = num./den;
tau = tau(:);

% calculate tangent t = r'/|r'|
tangent = [dx, dy, dz];
tangent = bsxfun(@rdivide, tangent,sqrt(sum(tangent.^2,2)));

% calculte binormal b = r' x r'' / |r' x r''|
binormal = [drxdrr1, drxdrr2, drxdrr3];
binormal = bsxfun(@rdivide, binormal,sqrt(sum(binormal.^2,2)));

% calcualte normal n = b x t
normal = cross(binormal, tangent);
normal = bsxfun(@rdivide, normal,sqrt(sum(normal.^2,2)));

end