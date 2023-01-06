function [Xout, Yout] = remove_initial_linear_cl(Xs,Ys)
% PURPOSE: This function addresses need to remove an initial section of the 
%          centerline. The model produces applies a constant upstream 
%          curvature boundary condition that results in a growing linear 
%          stretch at the mouth of the centerline. This linear stretch 
%          should be removed before computing further statistics.

% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com

% INPUTS: Xs, Ys      - coordinates of channel centerline

% OUTPUTS: Xout,Yout  - centerline coordinates with interpolated nodes


% threshold below which points are considered part of initial straight line 
% that is a product of the upstream boundary condtion
line_thresh = 1; 
z = 2; % how many points to use to fit the line
% find the slope and intercept of the line passing through z points
[fitline] = polyfit(Xs(1:z),Ys(1:z),1); 
% find difference between line and actual centerline
diffs = Ys - (fitline(2)+fitline(1)*Xs); 
% find centerline nodes that are not on the line
bigenough = abs(diffs) > line_thresh;  
firstone = find(bigenough==1,1); % identify the first node not on the line
Xout = Xs(firstone:end);
Yout = Ys(firstone:end);

end
