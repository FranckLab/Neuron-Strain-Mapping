function rf = findSeg(tree)
%
% findSeg(tree)
%
% Finds the indices of each segment that makes up the neuron. A segment is
% defined as either endpoint to soma or endpoint to branch point. 
% 
%
% Returns a cell size 1 x n, n is number of neurons. Each cell contains an
% array where each row represents a the indices of the nodes that make up a
% segment of the corresponding neuron. Array is padded with zeros for
% matrix agreement.
%
% INPUTS
%------------------
% tree:             Cell array of tree structures
%--------          
% OUTPUTS            
% -----------------
% rf:               Cell array of matrices containing the indices to nodes
%                   that make up segments of a neuron. A matrix with be an
%                   m x n matrix where m is the number of segments and the
%                   nonzero elements in n are the node indices. The matrix
%                   is padded with zeros for dimension agreeement. 
%

% Check tree is a cell
if ~iscell(tree), tree = {tree}; end

% Pre-allocate 
rf = cell(length(tree),1);

for jj = 1:length(tree) % number of neurons
    
    % Returns indices of branch points
    B = B_tree(tree{jj});
    Bx = find(B);

    % Gives path to parent for each node. (size(ipar,1)=number of points
    % that make up the neuron) (The last nonzero element in each row == 1,
    % which represents the soma)
    ipar = ipar_tree(tree{jj});

    % Gives indices of termination points
    T = T_tree(tree{jj});
    Tx = find(T);

   % Extracts the rows from ipar that begin with a termination point and
    % transposes it so that each column of df represents one full segment of
    % the neuron.
    df = ipar(Tx,:)';

    % Takes each column of df and makes it one really long column, but gets rid
    % of all the zero values elements that ipar_tree makes to keep the matrix a
    % valid size (number of elements vary so zeros are fillers)
    dfv = nonzeros(df);

    % Returns logical vector showing where the branch points are
    Lia = ismember(dfv,Bx);

    % Picks out the locations of the branch points. Reason for adding 1 to the
    % top of the vector is for indexing in the upcoming loop.
    ind = find(Lia);
    ind = [1;ind];

    % Pre-allocate. Each row of r represents a full segment (end to fork, end
    % to soma, or fork to soma). r is padded with zeros. 
    r = zeros(length(ind)-1,max(diff(ind)));
    
    for ii = 1:length(ind)-1 % number of segments

        % Determines the number of elements that make up each full segment
        q = (ind(ii+1) - ind(ii)) + 1;

        % Extracts a segment
        extract = dfv(ind(ii):ind(ii+1));

        % Sometimes the soma comes up at the beginning and end of a segment so
        % the value of 1 at the beginning is removed. 
        if extract(1) == 1
           q = q - 1;
           r(ii,1:q) = extract(2:end);
        else
           r(ii,1:q) = extract;
        end

    end

    % Removes repeated segments
    r = unique(r,'rows');
    
    % Remove rows with 3 or less nonzero elements because the third
    % derivative must be found later to calculator torsion
    check = r;
    check(check>0) = 1;
    check = sum(check,2);
    r(check<4,:) = [];

    % Sort each row in accending order
    r = sort(r,2);

    rf{jj} = r;

end % end number of neurons

end % end function

