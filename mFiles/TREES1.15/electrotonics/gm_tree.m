% GM_TREE   Membrane conductances of the segments of a tree.
% (trees package)
%
% gm = gm_tree (intree, options)
% ------------------------------
%
% returns the membrane conductance of all elements [in Siemens].
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - options::string: {DEFAULT: ''}
%     '-s'  : show
%
% Output
% -------
% gm::Nx1 vector: membrane conductance values of each segment
%
% Example
% -------
% gm_tree (sample_tree, '-s')
%
% See also gi_tree
% Uses surf_tree ver_tree Gm
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function gm = gm_tree (intree, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use only membrane conductance vector/value for this function
if ~isstruct (intree),
    Gm = trees {intree}.Gm;
else
    Gm = intree.Gm;
end

if (nargin < 2)||isempty(options),
    options = ''; % {DEFAULT: no option}
end;

gm = Gm .* surf_tree (intree) / 100000000;
% conversion from um2 to cm2

if strfind (options, '-s'), % show option
    ipart = find (gm ~= 0); % single out non-0-length segments
    clf; shine; hold on; plot_tree (intree, gm, [], ipart); colorbar;
    title  (['membrane conductances (total: ' num2str(sum (gm)) ...
        ') [S]' ]);
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view(2); grid on; axis image;
end

