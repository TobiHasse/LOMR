function [stream_length, tortuosity] = quicktor_TRH(Xs,Ys)
% PURPOSE: This function calculates the tortuosity in a two dimensional
%          line of many points

% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com
%          Tobias Hasse, March 2015 updated to improve speed

% INPUTS: Xs, Ys         - coordinates of channel centerline

% OUTPUTS: stream_length - cumulative length of stream centerline
%             tortuosity - the ratio of stream length to valley length

% remove linear section of stream at mouth imposed by boundary condition
[Xsf,Ysf] = remove_initial_linear_cl(Xs,Ys); 

stream_length = sum(sqrt(diff(Xsf).^2+diff(Ysf).^2));  %TRH new method

val_length = sqrt((Xsf(1)-Xsf(end))^2 + (Ysf(1) - Ysf(end))^2);
tortuosity = stream_length/val_length;

end
