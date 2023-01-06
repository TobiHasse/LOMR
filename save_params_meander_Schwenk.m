function save_params_meander_Schwenk()
% Purpose:  Provide commentary for the Schwenk input parameters and
%           describe adjustments to these parameters
%           Create this file in commented code for better communication and
%           understanding
%           Many of these variables are for the Beatton River, Canada as
%           published by Parker and Andrews (1985)
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June, 2021

fprintf('Welcome to: "%s" \n',...
        fullfile(fileparts(mfilename('fullpath')),mfilename))


B      = 35;        % half width of the channel (meters)
Qo     = 325.6;     % Characteristic streamflow (cubic meters / second) 
Eo     = 1.85*10^-8;% erodibility coefficient 
Cfo    = 0.0036;    % friction factor, affects meander bend length
alphaao= 10;        % A, the cross slope factor.  This value has 1 added to 
                    % it within the flowfield computation, this is an error
                    % in flowfield.  See Tobias Hasse's dissertation

% the following values for depth and velocity are not used, and are 
% overwritten in the meander model.  They correspond to a valley slope of: 
% while this error does not appear to have affected the model it could have
% affected other calculations that were made based on Schwenk's input file:
% params.mat.  Any potential errors did not affect Hasse's dissertation
S_valo = 0.0067;    % slope of the valley published by Schwenk (2015)
Do = 1.0582;        % bank full water depth
Uo = 4.3955;        % innitial reach averaged bank full water velocity
                    
S_valo = 0.00067;   % slope of the valley used by Schwenk and Hasse for 
                    % simulation
                    
% node spacing thresholds upper and average limits (meters)
dS_spacing_thresh  = 105;
avg_spacing_thresh = 1.2;

%% computing the correct and incorrect depth and velocity
g          = 9.81;  % gravity m/s^2

% incorrect depth and velocity
S_valo_incorrect = 0.0067               % INCORRECT valley slope
Do_incorrect     = (Qo/2/B)^(2/3)*(Cfo/g/S_valo_incorrect)^(1/3) 
Uo_incorrect     = Qo/(Do_incorrect*B*2)  

% correct depth and velocity
S_valo_correct = 0.00067                % CORRECT valley slope
Do_correct     = (Qo/2/B)^(2/3)*(Cfo/g/S_valo_correct)^(1/3) 
Uo_correct     = Qo/(Do_correct*B*2)     



save('params.mat','B','Qo','Do','Uo','S_valo','Eo','Cfo','alphaao',...
    'dS_spacing_thresh','avg_spacing_thresh','-v7.3')

end