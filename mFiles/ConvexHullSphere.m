function I = ConvexHullSphere(L1,R1,L2,R2,um2vxl)
%
% ConvexHullSphere(L1,R1,L2,R2,um2vxl)
%
% Create the convex hull of two spheres
%
% INPUTS
%------------------
% L1:           Center of first shere in non-dimensional units 
% R1:           Radius of first sphere in non-dimensional units
% L2:           Center of second shere in non-dimensional units 
% R2:           Radius of second sphere in non-dimensional units
% um2vxl:       1x3 vector specifying the x,y,z conversion factors
%--------
% OUTPUTS
% -----------------
% I:            Logical volume with convex hull mask of the two spheres
%

% Scale the z-step
scale = um2vxl(3)/um2vxl(1);

% Reorder the coordinates to y x z
L1 = [L1(2) L1(1) L1(3)];
L2 = [L2(2) L2(1) L2(3)];

% Find the length of the line between both points
d = norm(L1-L2);

% Unit vector 
nhat = (L2-L1)/d;

% Diameter of the spheres (run)
res = 10000;
x1 = linspace(-R1,R1,res); x2 = linspace(-R2+d,R2+d,res);

% Take the axisymmetric version (rise) 
y1 = sqrt(R1^2-x1.^2)/scale; y2 = sqrt(R2^2-(x2-d).^2)/scale;

% tan1 thus equals tan2 on the respective surfaces 
% and you get the boundary tangent for revolution
tan1 = -x1./y1;

% try connecting with the x pt with the same slope, minimizing the error bt.
% the line segment and the tangent slopes
segslp = (y2-y1)./(x2-x1);
[~, ind] = min(abs(segslp-tan1));
% tangent line equation is: y-y1(ind) = tan1(ind)*(x-x1(ind))

% find the height limits for "in the hull" off the colinear axis starting at
% one tangent pt and ending at the other
h1 = y1(ind); h2 = y2(ind);

% shrink the box to the spheres, leave everything else as a zero
Lim = [min(L1-R1,L2-R2); max(L1+R1,L2+R2)];
[mx, my, mz] = meshgrid(Lim(1):Lim(2),Lim(3):Lim(4),Lim(5):Lim(6)); 
m = [mx(:), my(:), mz(:)];

%return to the original setup, now with the hlim function that takes care of the
%"between the spheres" values
newijk = bsxfun(@minus,m,(L1 + x1(ind)*nhat));
dot_ = sum(bsxfun(@times, newijk, nhat),2);
distmx_ = newijk - bsxfun(@times, nhat,dot_);
distmx = sqrt(sum(distmx_.^2,2)); 
distmx = reshape(distmx,size(mx));
linefrac = dot_/(x2(ind) - x1(ind));
linefrac = reshape(linefrac,size(mx));

% Create a logical volume
ConvexHull = linefrac >= 0 & linefrac <= 1 & (h1+(h2-h1)*linefrac) >= distmx;
Sphere = sqrt((mx - L1(1)).^2 + (my - L1(2)).^2 + (mz - L1(3)).^2/scale^2) <= R1;
I = ConvexHull | Sphere;

end