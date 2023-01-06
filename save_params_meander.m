function save_params_meander()
% Purpose:  Organize many of the adjustable inputs for the meander
%           migration model as used by Tobias Hasse in his dissertation
%           Create this file in commented code for better communication and
%           understanding
%           Many of these variables are for the Beatton River, Canada as
%           published by Parker and Andrews (1985)
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June, 2021

B      = 35;        % half width of the channel (meters)
Qo     = 325.6;     % Characteristic streamflow (cubic meters / second) 
S_valo = 0.00067;   % slope of the valley
Eo     = 3.7*10^-8; % erodibility coefficient (double the value published 
                    % for the Beatton River to speed up migration rates)

% node spacing thresholds upper and lower limits (meters)
dS_spacing_thresh = 133;
too_close_thresh = 19.95;

save('params_meander.mat','B','Qo','S_valo','Eo','dS_spacing_thresh',...
    'too_close_thresh','-v7.3')

end