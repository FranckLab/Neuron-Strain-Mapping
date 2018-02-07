function tree = resampleTree(tree)
%
% resampleTree(tree)
%
% Make nodes within tree structure equidistant by a value of 25th
% percentile of all segment lengths making up the neuron
%
% INPUTS
%------------------
% tree:          tree structure. Can be either cell array or a single
%                structure.
%--------            
% OUTPUTS            
% -----------------
% tree:          Resampled tree structure as a cell array even if input was
%                originally a structure

% Check tree is a cell
if ~iscell(tree), tree = {tree}; end

for ii = 1:length(tree); % Number of neurons

    % Find the lengths of each segment in micrometers
    len = len_tree(tree{ii});

    % 25th pecentile of the len output
    Y = prctile(len,25);

    % Resampling distance is chosen as 25th percentile of all edge lengths
    tree{ii} = resample_tree(tree{ii},Y);

end

end

