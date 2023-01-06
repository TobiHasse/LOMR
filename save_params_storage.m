function save_params_storage(dt_save_years)
% Purpose:  organize many of the adjustable inputs to the storage time
%           distributins and put them in a parameters file which can be
%           loaded rather than having a basket full of variables to hand
%           between functions.
%           Create this file in commented code for better communication and
%           understanding
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June, 2021


% see starts_and_ends.m for better starting and ending points, the
% following were used for Tobias Hasse's dissertation
x_start = 345;      % x pixel coordinate for start of upstream 
                    % reach for computing storage time on floodplain
x_end   = 2240;     % x pixel coordinate for end of downstream 
                    % reach for computing storage time on floodplain
end_sim = 6845;     % ending simulation before final planform available 
                    % due to the river getting unusually short near the 
                    % end of the simulation

dt = 15;            % time step between saved river planforms which 
                    % are 2 years apart

time_step_years   = dt * dt_save_years;   
sim_time_ky       = round(end_sim * time_step_years / 1000);

% time step to display progress of stratigraphic model and channel
% occupation of floodplain pixels. 5000 years recommended
display_progress_years = 5000; 

pix_per_chan      = 5;    % width of channel in pixels

implicit_age_dist = true; % this sets the age distributions to be computed 
                          % based on the storage time distributions
save_images       = true; % saves intermediate images in display_model.m
load_limits       = false;% load boundaries of subsections as for 
                          % dissertation.  Must have right amount of RAM,
                          % see array_size below

% calculating storage times requires a _lot_ of RAM, using virtual memory
% slows things down significantly (on spinning hard drives) setting the
% array size for the sedimentary stack to be a portion of the available RAM
% can help minimize virtual memory use and speed up run times.  
% array_size also sets the size of subsection of the floodplain to analyze
% dynamically set array size:
hm = memory; % this works on a PC but might not on a MAC
array_size = hm.MaxPossibleArrayBytes/45; 
% manually set array size:
%     array_sz = 10*10^6; % PC 4 GB RAM\\ depends on machine memory
%     array_sz = 220*10^6; % PC 8 GB RAM\\ depends on machine memory
%     array_sz = 450*10^6; % Mac 16 GB RAM\\ this for dissertation
array_size = 100*10^6;  % optionally set array size                                                 

save('params_storage.mat','x_start','x_end','end_sim','dt',...
    'dt_save_years','time_step_years','sim_time_ky',...
    'display_progress_years','pix_per_chan',...
    'implicit_age_dist','save_images','load_limits','array_size','-v7.3')


end