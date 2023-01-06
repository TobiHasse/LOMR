function [Xs, Ys, cutidcs,t1,t2,t3] = TRH_rm_cutoffs(Xs, Ys,...
    search_radius, R)

% This function steps downstream along a centerline seeking non-adjacent
% centerline nodes within a specified search distance.  If any are found
% the node at a minimum distance is chosen and all the nodes between the
% two are removed.  Node removal is done in the upstream direction.

% This is adapted by Tobias Hasse April 2015 tobiack@udel.edu, from 
% Jon Shwenk remove_cutoffs.m, 8/20/2014. jonschwenk@gmail.com

% TRH_rm_cutoffs.m replaces the functions 
% cutoffs.m intersections.m and gen_banks.m which cause bugs in
% smooth_after_cutoffs.m
% TRH_rm_cutoffs.m is not dependant on the Statistics and Machine Learning
% Toolbox for MATLAB

% INPUTS: Xs, Ys - coordinates of pre-cut channel centerline
%         search_radius  - Maximum distance for 2 nodes to create a cutoff
%         R              - Range of nodes to skip when searching for cutoff

% OUTPUTS: Xs, Ys - channel centerline coordinates with cutoffs removed
%          cutidcs - Mx2 array, where M = number of cutoffs identified.
%                    cutidcs(M,1) = upstream cutoff node 
%                    cutidcs(M,2) = downstream cutoff node.

% t1=0; t2=0; t3=0;
% tic

cutidcs = [];
N = numel(Xs);

%% method that tests for cutoff every time
i=1; j=0;                
% while i < (N-R)
%     dist = ((Xs(i+R:end)-Xs(i)).^2+(Ys(i+R:end)-Ys(i)).^2).^.5;
%     [mn id_cut] = min(dist);
%     if mn < search_radius
%         cutidcs=[cutidcs; i   id_cut+i+R-1];  % cuts out closest nodes
% %       cutidcs=[cutidcs; i+1 id_cut+i+R-2];  % this cuts fewer nodes
%         j=j+1;
%         i = id_cut+i+R-1;             % jump downstream past cutoff nodes
%     end
%     i = i +1;
% end
% for i = numel(cutidcs)/2:-1:1
%     Xs(cutidcs(i,1):cutidcs(i,2))=[];
%     Ys(cutidcs(i,1):cutidcs(i,2))=[];
% end
%   

%% new method, only tries to cut if there are any cutoffs, 
% no testing in big loop
% this method is 70% longer (170%) compared with Shwenk method (presumably
% the method using RangeSearch from the Statistics & Machine Learning tools
% t1=toc; tic
mn     = zeros(N-R,1);
id_cut = zeros(N-R,1);
% instead of computing until end, create a ds_idx that is some long
% distance downstream
% or make vector of the starting index, then repmat and do matrix algebra?
for i = 1:(N-R)
%     distance = ((Xs(i+R:end)-Xs(i)).^2+(Ys(i+R:end)-Ys(i)).^2).^.5;
    [mn(i) id_cut(i)] = min(((Xs(i+R:N)-Xs(i)).^2+...
        (Ys(i+R:N)-Ys(i)).^2).^.5);
%     [fmn(i) fid_cut(i)] = find(((Xs(i+R:end)-Xs(i)).^2+...
%         (Ys(i+R:end)-Ys(i)).^2).^.5)<thresh;
%     [~,mn_id] = min(fmn);
%     id = fid_cut(mn_id);
end

i = 1;
search_id=find(mn<search_radius,1);  % find the most upstream cutoff
if search_id > 0                     % add cutoff idcs to cutidcs
    TRA_rm_cutoffs_slice_n_dice___L=64;
    cut_idx = search_id+id_cut(search_id)+R-1;
    cutidcs = [cutidcs; search_id cut_idx];
    i = i + cut_idx;
    % returns how many nodes farther downstream the next cutoff starts
    search_id=find(mn(i:end)<search_radius,1); 
    while search_id                  % add more cutoffs only if they exist
%         id_cut(i+search_id-1);
        cut_idx = search_id+id_cut(i+search_id-1)+i+R-2;
        cutidcs = [cutidcs; search_id+i-1 cut_idx];
        i = cut_idx;
        % returns how many nodes farther downstream the next cutoff starts
        search_id=find(mn(i:end)<search_radius,1); 
    end
    for i = numel(cutidcs)/2:-1:1
        Xs(cutidcs(i,1):cutidcs(i,2))=[];
        Ys(cutidcs(i,1):cutidcs(i,2))=[];
    end
    
end
% keyboard
% t2=toc;






