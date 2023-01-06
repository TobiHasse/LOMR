function [Xs, Ys, cutidcs, crofset] = remove_cutoffs_TRH(Xs, Ys, ...
    search_radius, R, crofset)
% This code works in a downstream direction along a set of centerline
% coordinates to detect and remove cutoffs by finding non-adjacent nodes 
% that are some threshold distance (search_radius) near each other.

% Jon Schwenk, 8/20/2014. jonschwenk@gmail.com
    
% Modified by Tobias Hasse April 2015, tobiack@udel.edu
% ************************************************************************
% See also discussion und "Spacing Criteria" in migration_model_TRH_Ch2.m
% ************************************************************************
% From fig 7b Schwenk et al 2015, the shortest cutoff Lcut/Lo = 10^1.2
% Lo appears to ~ 9B, so the streamlength of the shortest cutoff = 135B
% Rather than checking some distance from the search node to exclude, I use
% the spacing thresholds to determine the number of nodes to skip: crofset.
% Since rangesearch returns a cell array of indecies, I shift those
% indicies by the row index + crofset.  If any value remains above 0 then
% there is a cutoff.  I use a simple algorithm to find all the cutoffs
% because efficiency doesn't matter (cutoffs are rare).  Since rangesearch
% returns indecies in sorted order, find the first element of the offsets
% is the correct downstream cut index.  
% crofset is passed around so it doesn't need to be recreated, only
% expanded occationally
    
% possible improvement (slower?) use near_elem_TRH to determine a distance 
% threshold for the downstream index to start checking for cutoffs.  
% Recalculate every time
% pseudocode:
% crofset= num2cell(near_elem_TRH(s+distance_thresh,s));

% INPUTS: Xs, Ys         - coordinates of pre-cut channel centerline
%         search_radius  - Maximum distance for 2 nodes to create a cutoff
%         R              - Range of nodes to skip when searching for cutoff
%         crofset        - cell array of row offsets

% OUTPUTS: Xs, Ys  - channel centerline coordinates with cutoffs removed
%          cutidcs - Mx2 array, where M = number of cutoffs identified.
%                    cutidcs(M,1) = upstream cutoff node 
%                    cutidcs(M,2) = downstream cutoff node.
%          crofset - passed back and forth to avoid recreating every step

% ensure column vectors
Xs = Xs(:);
Ys = Ys(:); 
cutidcs  = [];       % must be returned, declare early even if un-needed

p = [Xs Ys];
[idx, ~] = rangesearch(p, p, search_radius,'Distance', 'euclidean', ...
    'NSMethod', 'kdtree'); % returns for each node all neighbors within 
                    % the search radius and their corresponding distances
N = numel(idx);

if N > numel(crofset)      % ensure offset vector is longer than id vector
    crofset=num2cell((R:(numel(idx)+R+500))');
    disp(['remove_cutoffs_TRH.m line: ',num2str(54),...
        ', extending row offsets to: ',num2str(numel(crofset)),' rows.'])
end

% a cutoff could occur IFF any offset_ids>0
offset_ids = cellfun(@minus,idx,crofset(1:numel(idx)),'UniformOutput',0);
cutoffs=max([offset_ids{:}]);

if cutoffs > 0 
% loop through and find the cutoffs only if they exist, since cutoffs are 
% rare, efficiency is not that important (± 0.1% of model run spent here)
    search_idx = 1; 
    % a while loop is used to skip downstream so there are no cutoffs from
    % nodes already cutoff
    while search_idx <= N
        id = find(offset_ids{search_idx}>0,1);   
        if id
            curr_row = idx{search_idx};
            cutidx = curr_row(id);
            cutidcs = [cutidcs; search_idx cutidx];
            search_idx = cutidx;             % jump downstream past cutoff
        end
        search_idx = search_idx + 1;
    end 

    % perform the cutoff (ie remove the nodes from the centerline)
    % Do this bottom up or indexing will be off if 2 cuttofs occur at the
    % same timestep
    for i = numel(cutidcs)/2:-1:1 
        Xs(cutidcs(i,1):cutidcs(i,2))=[];
        Ys(cutidcs(i,1):cutidcs(i,2))=[];
    end
end % if
end % function
