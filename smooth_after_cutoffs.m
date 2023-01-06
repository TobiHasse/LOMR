function [Xs, Ys, cutidcs] = smooth_after_cutoffs(Xs, Ys, cutidcs, Npre)
% PURPOSE: This function will smothe sharp bends in curvature which occur 
%          after cutoff.  

% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com

% INPUTS: Xs, Ys  - coordinates of channel centerline
%         cutidcs - the indices where cutoffs occured
%         Npre    - the number of points in the centerline before any 
%                   cutoffs were removed.
% OUTPUTS: Xs,Ys  - centerline coordinates smothed cutoffs
%         cutidcs - the indices where cutoffs occured

numsmooth = 15;%number of nodes before and after cutoff to use in smoothing
SGnfilt = 1; % number of times to filter interpolated coordinates
SGorder = 5; % order of the filter for interpolated coordinates -> lower 
%              order suppresses curvature more
add_pts = 2; %number of points to add to interpolated section across cutoff

% keyboard

ncuti = 0;
smoothidcs = cutidcs(:,1)-1;
%this loop finds the index of the node immediately upstream of each cutoff; 
% indices must be adjusted for upstream cutoffs
for i = 1:numel(smoothidcs) 
    if i == 1
    else
        ncuti = ncuti + (cutidcs(i-1,2)-cutidcs(i-1,1)+1);
        smoothidcs(i) = smoothidcs(i) - ncuti; 
    end
end 

M = numel(smoothidcs);
for cutno = M:-1:1 % work upstream to downstream
    % if the cutoff occurs at the beginning of the centerline
    if cutidcs(cutno,1) == 1 
        % if the cutoff occurs at the end of the centerline
    elseif cutidcs(cutno,2) == Npre 
    else        
               
    cutoff = smoothidcs(cutno);
    cuton = cutoff + 1;
    
    % interpolate two nodes linearly across the cutoff
    x1 = Xs(cutoff); x2 = Xs(cuton);
    y1 = Ys(cutoff); y2 = Ys(cuton);
    
    dx = x2 - x1;
    dy = y2 - y1;
    
    xi1 = x1 + (1/3)*dx;
    xi2 = x1 + (2/3)*dx;
    yi1 = y1 + (1/3)*dy;
    yi2 = y1 + (2/3)*dy;
    % in case the cutoff occurs near the beginning of the river
    numsmooth_beg = min(numsmooth, cutoff-1);   
    % in case the cutoff occurs near the end of the river
    numsmooth_end = min(numsmooth, length(Xs)-cuton-1); 
    % have to convert to double to suppress a warning about mixed data 
    % types from ode45 in interparc.
    Xsmooth = double([Xs(cutoff-numsmooth_beg:cutoff); xi1; xi2;...
        Xs(cuton:cuton+numsmooth_end)]); 
    Ysmooth = double([Ys(cutoff-numsmooth_beg:cutoff); yi1; yi2;...
        Ys(cuton:cuton+numsmooth_end)]);

    % Interpolate across the cutoff
    % s_pts is a function of the streamlength, see further in code. 
    % The initial ratio of pts/length is maintained.
%   [XY] = interparc(numsmooth*2+add_pts,Xsforsmooth,Ysforsmooth,'spline'); 
%     Xsmooth = single(XY(:,1)); % convert back to a single 
%     Ysmooth = single(XY(:,2));

    % Use Savitzky-Golay filter to smooth interpolated coordinates around 
    % cutoff window must be larger than order of filter, hence the max(). 
    % must also be odd, hence the rounding.
    window =  max(2.*round((SGorder+1+1)/2)-1,2.*round((length(Xsmooth)/...
        4+1)/2)-1); 
    tmp     =   [Xsmooth, Ysmooth];
    for iFilt=1:SGnfilt;
        tmp =   sgolayfilt(tmp,SGorder,window);
    end
    Xsmooth = tmp(:,1);
    Ysmooth = tmp(:,2);
    
    Xs = [Xs(1:cutoff-numsmooth-1); Xsmooth; Xs(cuton+numsmooth+1:end)];
    Ys = [Ys(1:cutoff-numsmooth-1); Ysmooth; Ys(cuton+numsmooth+1:end)];    
    end
    
    cutidcs(cutno,2) = cutidcs(cutno,2) + add_pts; %account for added nodes
    % account for cutoffs near the end but not exactly at end of centerline
    if cutidcs(cutno,2) > numel(Xs) 
        cutidcs(cutno,2) = numel(Xs);
    end
        
end % for

end % function

% clf
% plot(Xs,Ys,'.');axis equal; hold on
% plot(Xs(cutoff),Ys(cutoff),'ro')
