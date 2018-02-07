function tree = taubinSmoothing(ri,tree)
%
% taubinSmoothing(ri,tree)
%
% Make nodes within tree structure equidistant by a value of 25th
% percentile of all segment lengths making up the neuron
%
% INPUTS
%------------------
% tree:          Tree structure. Can be either cell array or a single
%                structure.
% ri:            Cell array containing the node indices of segments. Use
%                output from findSeg for this input. 
%--------            
% OUTPUTS            
% -----------------
% tree:          Smoothed tree structures stored in cell arrays. 


% Check ri and tree are cells
if ~iscell(ri), ri = {ri}; end
if ~iscell(tree), tree = {tree}; end

% Constants
lambda = 0.63;
mu = 1/(.1-(1/.63)); % According to equation

for kk = 1:length(tree);
    
    % Get segments from one neuron at a time
    r = ri{kk};

    for ii = 1:size(r,1)

        % Extract single segment
        a = r(ii,r(ii,:) > 0);

        % Find corresponding coordinates for each index
        c = [tree{kk}.X(a) tree{kk}.Y(a) tree{kk}.Z(a)];

        % Pre-allocate
        dlap = zeros(size(c,1)-2,3);

        for n = 1:100

            % Find discrete Laplacian
            for jj = 2:size(c,1)-1        
                dlap(jj-1,:) = 0.5*(c(jj-1,:) - c(jj,:)) + ...
                    0.5*(c(jj+1,:) - c(jj,:));    
            end

            if mod(n,2) == 0
                % Shrinking step
                c(2:end-1,:) = c(2:end-1,:) + lambda*dlap;
            else
                % Expansion step
                c(2:end-1,:) = c(2:end-1,:) + mu*dlap;
            end    

        end

        % Save smoothed values back into trees
        tree{kk}.X(a) = c(:,1);
        tree{kk}.Y(a) = c(:,2);
        tree{kk}.Z(a) = c(:,3);

    end % end smoothing

end % end number of neurons

end % end function

