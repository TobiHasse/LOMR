function [angles2] = cl_angles(Xs, Ys)
% PURPOSE: calculate the angles of each centerline node relative to the
%          donstream direction of flow.
% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com

% INPUTS: Xs, Ys    - coordinates of channel centerline

% OUTPUTS: angles2  - angle at each node

% Xs=X;Ys=Y;
% Reshape inputs
Xs = Xs(:)';
Ys = Ys(:)';
% if size(Xs,1)~=1    
%     Xs = Xs';
% end
% if size(Ys,1)~=1
%     Ys=Ys';
% end   

% Estimate valley centerline via linear trendline through channel
% centerline. 
% v_cl = polyfit(Xs,Ys,1); % fit a linear trendline as initial valley 
% direction guess

% See http://www.mathworks.com/matlabcentral/newsreader/view_thread/276582
% for method of calculating angle between two lines/vectors. Angle
% calcuated in radians, range is from -pi to pi.
% for u = 1:length(Xs)-1
%     dxx = Xs(u+1)-Xs(u);
%     V1 = [Xs(u)+dxx, Ys(u)+(v_cl(1)*dxx)];
%     V2 = [Xs(u+1)-Xs(u), Ys(u+1)-Ys(u)];
%     if (V2(1)<0 && V2(1)>0) % quadrant 2
%         angles(u) = atan2((Ys(u+1)-Ys(u)),Xs(u+1)-Xs(u))+pi; % IV
%     elseif V2(1)<0 && V2(1)<0 % quadrant 3
%         angles(u) = atan2((Ys(u+1)-Ys(u)),Xs(u+1)-Xs(u))+pi; % IV
%     else
%         angles(u) = atan2((Ys(u+1)-Ys(u)),Xs(u+1)-Xs(u));
%     end
%         angles2(u) = atan2((Ys(u+1)-Ys(u)),Xs(u+1)-Xs(u));
% end
Xang = Xs(1:end-1);
Xang_plus1 = Xs(2:end);
Yang = Ys(1:end-1);
Yang_plus1 = Ys(2:end);
angles2 = atan2(Yang_plus1-Yang,Xang_plus1-Xang);

% Angle at first point cannot be estimated from coordinates, so backwards-
% fit a spline to interpolate it. assumes some kind of continuity of angles
% along the streamwise direction

% Linear interpolation for first angle
ang_first = angles2(1) - (angles2(2)-angles2(1));
angles2 = [ang_first angles2]';

% Spline interpolation for first angle
% x_first = [2:length(angles)+1];
% ang_first = spline(x_first, angles, 1); 
end