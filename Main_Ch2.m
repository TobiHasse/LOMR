%% Readme (and read the supplementary Readme if you plan on modifying code)
% This script runs three components as described in the accompanying paper:
% a meander migration model, an algorithm to extract individual bends, and 
% finally an algorithm to plot the results.

% The meander migration model is run for the parameters given in the paper,
% but with a shorter initial centerline. You can lengthen it if you'd like,
% but I wanted to keep runtime on the order of minutes instead of hours for
% this example.

% I hope you get under the hood of the model and extraction algorithm, but
% if you play with parameters or node spacings, it is likely that you will
% encounter numerical instability. Instability will often cause problems
% in the 'enforce_spacing' subroutine, and I have added a warning and
% 'keyboard' command for debugging if that happens. 
%% generate parameter input files
% rather than putting parameters directly in *.mat files, I have created
% commented code which will generate input parameter files.  
% This allows you to see commentary and variable definitions, or
% make adjustements for your simulations

close all
clear
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% cd 'C:\Users\thasse\Documents\MATLAB\meander code\Ch2\bug fix' %
cd 'C:\Users\thasse\Documents\MATLAB\test' %

dt_save_years = 2;                  % save river planform every ## years

save_params_meander()               % parameters for meander model
% save_params_meander_Schwenk()     % as used by Schwenk (with commentary)

save_params_storage(dt_save_years)  % for computing storage times
save_params_deposition()            % sediment deposition 
save_initial_planform()             % initial random planform 
                                    % included for (reproducibility)

%% initialize variables for meander model

do_you_have_stats_toolbox = 0;      % enter 1 if you have the Statistics 
                                    % Toolbox, else enter 0
sim_time_ky = 15; %211;              % simulation time in thousands of years
                                    % note orig sim deleted after 205.55
sim_time = sim_time_ky * 1000;
dt_save = 12;                       % save every ## iterations
dt =     dt_save_years / dt_save;

%% Run the meander migration model


% Edit this code block ***********************************************

% outfile = sprintf('Ackerman_3chan_%ska',num2str(sim_time_ky))
% outfile = sprintf('Hasse_3chan_120-130 %ska',num2str(sim_time_ky))
% outfile = 'double check 2016 params'
% [river, B, mxub, freq_mig] = check2016run( do_you_have_stats_toolbox,...
%     outfile,sim_time,dt,dt_save); %

% i=1; CFO = 0.0036; Ain = 10;
% outfile = 'Schwenk low slope'
    
% i=1; CFO = 0.024; Ain = 3;
% outfile = 'Hasse 2016 low slope'
% [river, nodecount, B, mxub, freq_mig] = migration_model_TRH(CFO(i),...
%     Ain(i), do_you_have_stats_toolbox,outfile,sim_time,dt,dt_save); %

% i=1; CFO = 0.024; Ain = 3;      % Hasse
% i=1; CFO = 0.01; Ain = 16;     % Testing
% i=1; CFO = 0.0036; Ain = 10; % Schwenk
i=1; CFO = 0.005; Ain = 12; % new → lambda = 11.5 ± 7.9

outfile = '15ka_2023' %Ch2 2016 params'
[river, nodecount, B, freq_mig] = migration_model_TRH_Ch2(CFO(i), ...
    Ain(i), do_you_have_stats_toolbox, outfile, sim_time, dt, dt_save); 

% outfile = 'double check 2016 params'
% [river, B, mxub, freq_mig] = check2016run( do_you_have_stats_toolbox,...
%     outfile,sim_time,dt,dt_save); %

%% Save output
% original save code:
% save(sprintf('%s_2yr_A3_Cfo24_2Eo.mat',outfile),'river','B','-v7.3') 
% more flexible save code:
save(sprintf('%s_2yr_A%d_Cfo%0.3f_2Eo.mat',...
    outfile,Ain,CFO),'river','B','-v7.3') 

%% Schwenk's visualization of output (optional)
% this script plays animations of the entire simulation and each cutoff
% bend
visualize_Schwenk


%% save output smaller (dissertation settings)
% take 2 year river, convert to 30 year 'riv' and 'riv2' (offset)
% dt = 15 for Hasse dissertation.  With 'river' saved every 2 years, dt = 
% 15 yeilds a 30 year time step for deposition
% load from file in case user changes dt for storage time computation
load params_storage dt 

start_time = 100 + dt ;
end_time   = numel(river); %50001; %25001
riv = river ([start_time:dt:end_time]);
% original save code
% save(sprintf('%s_30yr_A3_Cfo24_2Eo.mat',outfile),'riv','B',...
%     'start_time','-v7.3') 
% more flexible save code
save(sprintf('%s_30yr_A%d_Cfo%0.3f_2Eo.mat',...
    outfile,Ain,CFO),'riv','B','start_time','-v7.3') 

% save output with 
start_time = 100 + dt+round(dt/2); 
end_time   = numel(river); %50001; %25001
riv2 = river ([start_time:dt:end_time]);
% original save code
% save(sprintf('Ackerman_211ka_30yr_A3_Cfo24_2Eo_offset.mat'),...
%     'riv2','B','start_time','-v7.3') 
% more flexible save code
save(sprintf('%s_30yr_A%d_Cfo%0.3f_2Eo_offset.mat',...
    outfile,Ain,CFO),'riv2','B','start_time','-v7.3') 



%% The following code is for analysing
% Meandering River Dynamics and Storage Time github repo: MRDAST
% MRDAST depends on the LOMR repo to generate the meandering river
% planform evolution and history





%% Read in the saved channel planforms and view

% ********** LOAD riv, B, riv2, TO RUN CODE BLOCKS BELOW ************

clear
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat')
% load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat')
% load('code cleanup_30yr_A16_Cfo0.010_2Eo.mat')
% load('code cleanup_30yr_A16_Cfo0.010_2Eo_offset.mat')
% load('code cleanup_30yr_A10_Cfo0.004_2Eo.mat')
% load('code cleanup_30yr_A10_Cfo0.004_2Eo_offset.mat')

% load('5ka_2023_30yr_A3_Cfo0.024_2Eo.mat')
% load('5ka_2023_30yr_A3_Cfo0.024_2Eo_offset.mat')
load('50ka_2023_30yr_A12_Cfo0.005_2Eo.mat')
load('50ka_2023_30yr_A12_Cfo0.005_2Eo_offset.mat')

beep; pause(1); beep



%% start and ending nodes of simulation

% estimated run time ::::::::::::::::::: 2 minutes

load params_storage.mat
[ x_start_rx, x_start_buffered, x_end_rx, end_sim_rx ] = ...
                            starts_and_ends( riv, riv2, pix_per_chan, B);

fprintf(['Farthest downstream starting pixel is %d but the suggested\n',...
    'upstream starting pixel is %d to buffer some odd meander bends',...
    ' near the starting node. \nStarting pixel for dissertation was',...
    '%d \n \n'], x_start_rx, x_start_buffered ,x_start)

fprintf(['Suggested end simulation model step is %d. \nEnd simulation ',...
    'for dissertation was %d\nWithin that time the farthest upstream',...
    ' pixel for the end of the channel was %d. \nLast pixel for',...
    ' dissertation was %d \n \n'], end_sim_rx, end_sim, x_end_rx, x_end )

%% Create figures and files showing an overview of the model 

% estimated run time ::::::::::::::::::: 18 hours

% This will create figures showing the topography, age, selected 
% stratigraphy, floodplain width, migration rate, etc. throughout the
% simulation.  It will also save some files including the figures, data
% files, and a data file map of the floodplain area which was visited by
% the channel at least once.

display_model(riv, riv2, B)

%% animate channel occupation 
% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

outfile = 'junk.gif'
slice_interval = 1;


[valley_w_pxl,occ_bool] = chan_occ(outfile,slice_interval);


%% slicing the model to manage memory demands

% estimated run time ::::::::::::::::::: 2 seconds
% Dependancy:   Must run after display_model.m has successfully completed
%               preferably there is a floodplain data file titled: 
%               "Previously Occupied by Channel. Time step 6845.mat"
%               you may need to edit the slices.m function or file name

% for a model domain containing about 2,000 pixels across the valley and
% 3,000 pixels down valley and 7,000 time steps the array containing that
% is 42 billion pixels, at double precision this is 40 to 80 GB of RAM !!!
% (per array). To manage these memory demands the model is sliced into 
% short reaches and sediment is built and the storage times measured.

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

load params_storage.mat
end_sim = end_sim_rx; 
show_slice_figures = true;

% show the slices
[x_starts, x_ends, ymn, ymx] = slices( array_size, end_sim, ...
    x_start, x_end, show_slice_figures);
%% Create deposits and calculate storage times

% estimated run time ::::::::::::::::::: 40 to 100     ****** DAYS *******

% cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other

restart = false
slice = 0; % zero for first slice

show_figures = false;           % figures require significant RAM and time

deposit_storage_time(riv, riv2, B, show_figures,  restart,...
     [], [], [], [], slice, [],[],[],[])


%% Restart deposit_storage_time after some slices have been calculated

% estimated run time ::::::::::::::::::: 40 to 100     ****** DAYS *******
% run time is the remainder of the run started above

% this functionality is useful if there is a power outage or crash or if 
% the program must be stopped and restarted.  This is set up to restart 
% after the most recently completed floodplain slice and the input file
% name must be updated. This program will ignore an incomplete slice.


clear
restart = true
show_figures = false;           % figures require significant RAM and time
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% load planform files
load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat')
load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat')
% load storage time files for the MOST RECENT slice
load distr_3chan_205k_30yr_3.mat    % "_##.mat"  UPDATE THE SLICE NUMBER ##

deposit_storage_time(riv, riv2, B, show_figures, restart,...
    x_starts, x_ends, ymn, ymx, slice, ...
    pt_bar_dist, vert_dist, vd_age_dist, pb_age_dist )


%% age distribution of area nearby to the channel

% estimated run time ::::::::::::::::::: 12-24 hours

lam_vc = 10.5 % characteristic meander wavelength (down valley axis)
mb_age_lim = 333; %333 for 9990 years ? 10,000 years
mb_start_step = 4000; % 4000 to start measuring the meander belt at 120kyr

meander_belt_age( lam_vc, mb_age_lim, mb_start_step, riv, riv2, B)


%% This script should be continued to call other scripts and functions 
% which include the analyses and create the figures for Chapter 2 of Tobias
% Hasse's dissertation

