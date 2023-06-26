function [river, nodecount, B, freq_mig ] = migration_model_TRH_Ch2(Cfo,...
    alphaao, statstoolbox,outfilename,sim_time_yrs,dt_yrs,save_dt)
%,Ub_hist_max,X,Y,inplanformname)
disp('welcome to Schwenks model modefied by Hasse')
% This function was written by Jon Schwenk 2014-2015 
% and is available here 
% http://onlinelibrary.wiley.com/doi/10.1002/2014JF003252/full
% along with other functions as supplemetary info for the article:
% Schwenk, J., Lanzoni, S., & Foufoula?Georgiou, E. (2015). The life of a 
% meander bend: Connecting shape and dynamics via analysis of a numerical 
% model. Journal of Geophysical Research: Earth Surface.
% This is the main function for computing the migration of a meandering
% channel
% Updates by Tobias Hasse, tobiack@udel.edu, April 2015 include:
% 1) calling updated functions that are further optimized (for speed &/or 
%    memory) from the Schwenk functions, 
% 2) commenting out unused variables using "%TRH" to differentiate from 
%    Schwenk comments
% 3) adding some variables to help the new functions be more efficient (
%    e.g.) cell_row_offset
% 4) and iserting timers used to time different sections of the loop.  
% Updated functions have the suffix _TRH and are available from Tobias 
% Hasse

% Note that the save_nodecount C section has been changed so that duplicate
% indecies are not saved.  This requires the updated atom_tracking_TRH
% function to be called from Model_Extract_and_Plot.

scaleup = 5.7;
clock_time_start = clock;
%% Load and format coordinates, parameters
load params_meander.mat; % load reachwide variables, parameters
load initialplanform % dt = 0.05 looked good, testing with dt = 0.2 (too 
%                      high), trying 0.1
% You can overwrite these parameters which are in the params*.mat file
% Cfo = .024;  % TRH update Beatton River CA data based on criteria k ~= 
%                epsilon o from Johanesson & Parker 1989 
% alphaao = 3;  % TRH ***** WARNING ***** alphaao has 1 added within 
%                 flowfield for unknown reasons (ref: Johannesson & Parker 
%                 1985). TRH: Schwenk used 10 not 3

% Eo = Eo*2;   % TRH the new Cfo, A values, slow down the migration rate
% **** The updated Eo is in the params_meander.mat file ****
X = scaleup*Xo; %(1:ceil(numel(Xo)/scaleup)); % to truncate scaled version 
                                              % to same domain space
Y = Yo;         %(1:ceil(numel(Xo)/scaleup));
% Generate initial centerline coordinates
% X = [0:0.5:150]'*B;
% Y = randn(numel(X),1);
x_max_orig = X(end);
% load par_x_y   % for loading X Y coordinates saved at a previous partial 
% save for testing based on a well developed planform 120kplanform is the 
% start of storage time analysis of the 2016 model run
% load 120kplanform 
% load 50kplanform    %to pick up where we left before...
% X=x50k;
% Y=y50k;
% clear x50k y50k
%% Boolean speed considerations: parallel and saving--what and how often?
% par_flow = 1;           % the flowfield can be parallelized
par_cut = 0;            % the intersections search can be parallelized 
                        % required by cutoff.m (no stats toolbox)
%TRH save_stats = 0;    % saves various statistics in structure called stat
save_riv_vars = 0;      % saves river variables Xcl, Ycl, mig_Xcl, mig_Ycl 
                        % (centerlines and migrations at each time step)
save_nodecount = 0;     % saves variables needed to perform node tracking 
                        % - used for atom identification
nodecount = 0;          %if not saving nodecount, gives something to return
% save_dt = 20;     %(TRH passed into function) to save every nth time step
% TRH approximate convolution integral by truncating its computation after 
% some threshold distance (can be hard-coded suggest 100B)
appx_conv_int = 1       
% TRH record distributions of ub, migration rates, etc
meas_ub_migration = 0;

save_mem = 0;        % uses single precision instead of double if activated
% save points number of points throughout simulation to save data (adds 
% time but prevents loss of data if the run fails for some reason)
%TRH save_points = 1; 
%TRH save_number = 0;       % used for accounting; do not change this 
                            % variable from 0.

%% Simulation variables  % timing commented out bc read into function
% sim_time_yrs  = 100000;   %(TRH passed into function) length of 
                            % simulation, years
% dt_yrs = 0.1;             % time step, years, larger dt causes problems 
                            % in enforce_spacing.m
sim_time = ceil(sim_time_yrs/dt_yrs); % number of time steps, no dimensions 
dt = dt_yrs*365.25*24*3600; % time step in seconds
% disp_progress =  10000;  % display t every disp_progress time step
disp_progress = floor( sim_time / 60 ); % show progress every 1/60th
% disp_progress = floor( sim_time / 2 );    % only show 50% complete model

% number of segments to divide channel width into; effectively sets 
% resolution of across-stream grid (minimum of 2; results unchanged for >2)
%TRH n_pts = 2;             

%TRH end_save_times = ceil(sim_time/save_points):ceil(sim_time/...
%                             save_points):sim_time;
%TRH end_save_times(end) = sim_time;
%TRH start_save_times = [0 end_save_times(1:end-1)]+1;

%% Prescribe boundary, initial conditions, constants, and parameters
% Constants
rho = 998.1;    % density of water, kg/m^3
g = 9.81;       % gravity, m/s^2

% Solve initial flow variables and parameters
[~, tortuosity] = quicktor_TRH(X,Y);% tortuosity calculated assuming 
                                    % straight valley
Cf_ch = Cfo;                        % initial friction factor
S_val = S_valo;                     % initial valley slope
W = B*2;                         % channel width, m, doesn't change in time
So = S_valo/tortuosity(1);          % initial stream slope
alphaa = alphaao - 1;               % alphaa doesn't change in time
% Update 6/24/23, previously alphaa was passed to flowfield*.m and modified
% the -1 modification is correct based on Johannessen & Parker 1985
% eqs 14 & 15 which include an ERROR resulting in some authors using +1

% initial reach-averaged depth (TRH is this only truely Do because 
% tortuosity ? 1?)
Do = (Qo/2/B)^(2/3)*(Cfo/g/So)^(1/3); 
Uo = Qo/(Do*W);                     % initial reach-averaged velocity
% half-width:depth ratio for straight channel characterized by valley slope
%TRH betaao = B/Do;                 % TRH (apparently unused)     
% initial (Froude #)^2 for straight channel characterized by valley slope
F2o = (Uo/sqrt(g*Do))^2;            
%% Spacing criteria - adjusting these could cause instabilities, other bugs

% ****** Given Chapter 3 of my (Tobias Hasse's) dissertation there is
% better theory to manage all of the spacing thresholds ********

% TRH search_excl_range speeds up remove_cutoffs_TRH.m
% search_excl_range is a number of nodes greater than search_radius
% donwstream, but less than the smallest cutoff from Schwenk 2015, fig 7
% TRH Here we have a bit of a rabbit hole:  
% Schwenk 2015, fig 7 is nondimensionalized using the length scale Lo
% Lo = Do/2Cfo, for Schwenk 2015 Cfo = 0.0036 and 
% Do = 2.282 using the equation above (1.0582 in his input params.mat file
% This means that Lo = 317, but Schwenk used 147 (personal comm Jan 2020)
% Using Lo = 317, and Schwenk figure 7, minimum Lcut is 150 B which is why
% I made the max search threshold much lower than this (at 50) to be 
% conservative
% Using Lo = 147, minimum Lcut is 70B, however, based on the
% adjustments made by Tobias Hasse (table 1 of dissertation) to create
% smaller bends minimum Lcut might be 46. 
% Fortunately the search exlcusion range was 22.8 which is still less than
% 46, but not as conservative.

% number of channel half-widths between centerline to detect cutoffs; 
% default value (2) is width of channel
search_radius = B * 2; 
search_excl_range = ceil(search_radius/too_close_thresh +2);
cell_row_offset = num2cell([search_excl_range:(numel(X)*4+...
    search_excl_range)]'); 
    % make the offsets vector big enough for some increse in nodes
if search_excl_range*dS_spacing_thresh/B > 30;  %50
    sprintf(strcat('remove_cutoffs.m is skipping %d channel half',... 
        'widths downstream, cutoffs can be as small as 46 half widths'),...
        ceil(search_excl_range*dS_spacing_thresh/B)) 
end
% TRH this approximates conv_int_thresh = 100 in flowfield_TRH
conv_node_thresh = ceil( 100 / dS_spacing_thresh * B * 7/4);  

%% Curvature boundary conditions, smoothing
bc1o = 0;           % upstream initial and boundary condition for curvature
bcset = 1;          % 1 for periodic, 2 for C(1)=0, 3 for small random 
                    % Gaussian, 4 for linear interpolation
num_C_smooth = 1;   % number of times to smooth curvature each iteration

%% Initialize model variables
D_ch  = nan(1,sim_time);    D_ch(1) = Do;
S_ch  = nan(1,sim_time);    S_ch(1) = S_valo;
U_ch  = nan(1,sim_time);    U_ch(1) = Uo;
F2_ch = nan(1,sim_time);    F2_ch(1) = F2o; 
%TRH optoinal for measuring reach wide parameters at each step
chan_len = nan(1,sim_time); bends = nan(1,sim_time); 
mig_mag_max = 3;
freq_mig = zeros( 1,31 );

%TRH a_p_2 = nan(1,sim_time);
%TRH a_p_4 = nan(1,sim_time);
%TRH a_p_6 = nan(1,sim_time);
bc1   = nan(1,sim_time);   bc1(1) = bc1o;
%TRH betaa = nan(1,sim_time); betaa(1) = betaao;

%% Implement single precision - quicker runs, smaller output files
% TRH this might not matter on 64 bit machines.  I think 64 bit machines 
% use the same amount of RAM for single and double precision
if save_mem == 1
    X = single(X);
    Y = single(Y);
    U_ch = single(U_ch);
    D_ch = single(D_ch);
    Cf_ch = single(Cf_ch);
    S_val = single(S_val);
    Qo = single(Qo);
    B = single(B);
    alphaa = single(alphaa);
    bc1(1) = single(bc1(1));
end
    
%% Initialize saved variables' structure arrays
if save_nodecount
    save_dt = 1;
    nodecount=[]; % make sure nodecount is blank before reassigning
    nodecount( ceil( sim_time / save_dt ) ).A_cutoff_rem = [];
    nodecount( ceil( sim_time / save_dt ) ).B_spacing_ins = [];
    nodecount( ceil( sim_time / save_dt ) ).C_duplicate_rem = [];
%TRH     cut_idx = cell(1,sim_time);
end

if save_riv_vars
    river( ceil( sim_time / save_dt ) ).Xmig = [];
    river( ceil( sim_time / save_dt ) ).Ymig = [];
end
river( ceil( sim_time / save_dt ) ).Xcl = [];
river( ceil( sim_time / save_dt ) ).Ycl = [];
%TRH optional for measuring wavelength of each bend only saved time steps
waves(numel(river)).length=[]; 

%% Begin modelling
% Initialize centerline  
% TRH if not saving riv_vars (Line) 179 then these are not preallocated
save_t = 1;
river(save_t).Xcl(:,1) = X;
river(save_t).Ycl(:,1) = Y;
% % timers for speed testing
% clock_curv = 0; clock_flow = 0; clock_angle = 0; clock_cutoff = 0; 
% clock_smooth = 0; clock_space = 0; clock_update = 0; tic
%%
algorithm_start = clock 

for t = 1:sim_time   % ****************************************************
% line_221_t_is = t 
% keyboard % continue with dbcont, quit with dbquit
tic
% TRH these statements appear to be redundant, and cause problems if saving
% only every nth time step.
% X = river(t).Xcl; % get centerline coordinates computed at 
% Y = river(t).Ycl; % previous time step

%% Solve hydrodynamics
% (dimensionalized) curvature and related variables
[C, S_cum, dS] = curvars_TRH(X,Y,bc1(t),bcset, num_C_smooth); 
% if bcset == 1, the upstream curvature at the next time step is equal to 
% the downstream-most calculable (non-interpolated) curvature value, 
% a la a periodic boundary condition
bc1(t+1) = C(end-1); 
% clock_curv = clock_curv + toc; tic
chan_len(t) = S_cum(end); 
% Schwenk 2014 flowfield solution
% [ub,Lo_in_Bs,tf1,tf2] = flowfield_Schwenk(B,C,D_ch(t),S_cum,alphaa,...
%     F2_ch(t),Cf_ch,appx_conv_int);

% Hasse 2015 flowfield_TRH solution (based on Schwenk's function)
% TRH, Schwenk's method using sinuous channel parameters for simulation 
% This method was used in Chapter 2 of Tobias Hasse's dissertation
% [ub] = flowfield_TRH(B,C,D_ch(t),S_cum,dS,alphaa,F2_ch(t),...
%     Cf_ch,appx_conv_int,conv_node_thresh); 

% TRH, correct method using hypothetical straightened reach parameters
% **** USE THE [ub] FUNCTION CALL BELOW, NOT ABOVE FOR CORRECT RESULTS ****
[ub] = flowfield_TRH(B,C,Do     ,S_cum,dS,alphaa,F2o     ,Cf_ch,...
    appx_conv_int,conv_node_thresh); 
% the above closeley approximates Parker & Andrews 1985 Beatton River

% TRH in the past I measured meander wavelength using fft on the curvature
% series, concatenating many streamlines together to make a long series the
% method below measures bends without fft, but rather uses ub
% ub is a smothed version of curvature, some small bends may be included
% bends is tracking the number of meanders at each t, this will be
% converted to meander wavelength later in this function
temp=zeros(size(ub)); 
temp(ub>0)=1; 
bends(t) = sum(abs(diff(temp)))/2; %meander bends, (pair: 1 left, 1 right)

% clock_flow = clock_flow + toc; tic %DEBUG Lo_mn = min(Lo_mn,Lo_in_Bs); 
% % uncomment these lines if not running both flowfields
% Ro = min(abs(1./C)); % minimum radius of curvature along reach
% nu0 = B/Ro;
% ub = flowfield_Zolezzi_u(X, Y, C, S_cum, dS, B, D_ch(t), Cf_ch(t),...
%     S_ch(t), nu0);
%% Perform migration given flow field

Eps_cl  = Eo*ub; % migration rate along the centerline

% Compute erosion along outer banks (n=1)
% angles computed under assumption of constant valley direction found via 
% linear regression
angle = cl_angles(X,Y);  
% TRH multiplying by B below is in error and is what was done by Schwenk
% and by Hasse in dissertation Chapter 2. 
% These lines remain for reproducibility of results
% ********** THE RESULT OF THIS ERROR IS THAT ACTUAL Eo IS B*Eo **********
dX_cl = -Eps_cl.*sin(angle)*U_ch(t)*dt*B; % dimensionalized dx
dY_cl = Eps_cl.*cos(angle)*U_ch(t)*dt*B; % dimensionalized dy
% ******** USE THE TWO LINES BELOW, NOT ABOVE FOR CORRECT RESULTS ********
% dX_cl = -Eps_cl.*sin(angle)*U_ch(t)*dt; % dimensionalized dx
% dY_cl = Eps_cl.*cos(angle)*U_ch(t)*dt; % dimensionalized dy

if meas_ub_migration %TRH measure migration rate of nodes
    % this appears to be much greater than the migration rate orthogonal 
    % to streamflow (TRH June 2016)
    mig_mag = sqrt( dX_cl.^2 + dY_cl.^2 );     
    figure(100)
        mh = histogram([mig_mag',mig_mag_max],length(freq_mig));
        freq_mig = freq_mig + mh.Values;
end

Xf = X; Yf=Y; %TRH for saving final X,Y coordinates that match ub
mig_Xcl = X + dX_cl; % migrated X coordinates
mig_Ycl = Y + dY_cl; % migrated Y coordinates
% clock_angle = clock_angle + toc; tic
%% Locate and perform cutoffs
Npre = numel(mig_Xcl); 
if statstoolbox
%     [X, Y, cut_idx] = remove_cutoffs(mig_Xcl, mig_Ycl, cutoff_thresh, B);
    [X, Y, cut_idx,cell_row_offset] = remove_cutoffs_TRH(mig_Xcl,...
        mig_Ycl, search_radius, search_excl_range, cell_row_offset);
else
%     fprintf(['WARNING: The function cutoff smooths cutoffs inside\nbut',...
%         ' using different parameters than smooth_after_cutoffs below\n'])
%     pause(2)
    % 2023 cutoff is cutting off when it should not. removing large
    % sections of river centerline
%     [X,Y,n_cutoffs,Xl,Yl,Xr,Yr,cut_idx,ncut] = cutoff(mig_Xcl,mig_Ycl,...
%         B,par_cut);
    % make TRH_rm_cutoffs call here for no stats pack, pbbly faster than 
    % cutoff b.c. similar speed to remove_cutoffs
    [X, Y, cut_idx] = TRH_rm_cutoffs(mig_Xcl, mig_Ycl, search_radius,...
        search_excl_range);
end
% keyboard
% smooth around cutoffs
% clock_cutoff = clock_cutoff + toc; tic
if isempty(cut_idx) == 0  % if there were cutoffs, smooth them
    [X, Y, ~] = smooth_after_cutoffs(X, Y, cut_idx, Npre);
    if save_nodecount
        nodecount(t).A_cutoff_rem = cut_idx;
    end
end    
% clock_smooth = clock_smooth + toc; tic
%% Enforce maximum streamwise spacing threshold 
% TRH combine enforce spacing [X,Y,interp_idcs,cut_idcs] = 
%                     enforce_spc(X,Y,dS_max_spc,dS_min_spc,save_nodecount)
[X, Y, interp_idcs] = enforce_spacing(X, Y, dS_spacing_thresh, ...
    save_nodecount);
if save_nodecount
    nodecount(t).B_spacing_ins = [interp_idcs]; 
end
%% Enforce minimum streamise spacing threshold
dS_ibn = sqrt(diff(X).^2+diff(Y).^2); 
% find indices where centerline nodes are too close
remove_duplicates = find(dS_ibn < too_close_thresh); 
if save_nodecount
    if isempty(remove_duplicates) == 0
        %TRH repmat is not needed for the updated atom_tracking, which
        %doesn't need two copies of the removed node at lines 132 & 135
        nodecount(t).C_duplicate_rem =  remove_duplicates; 
        %TRH sort(repmat(remove_duplicates,1,2)); 
        % the sort(repmat()) is because each adjustment of nodes has a 
        % starting and ending node, so if it's only one node being removed 
        % there should be two identical entries in the structure
    end
end
X(remove_duplicates)=[];
Y(remove_duplicates)=[];
% clock_space = clock_space + toc;tic
%% Store variables
if ~mod(t,save_dt) %TRH saves river planform, but not every time step
    save_t = save_t + 1;  % increment river save counter
    river(save_t).Xcl = X;
    river(save_t).Ycl = Y;
    % TRH lenngth of each bend to left or right ÷ half width and tortuosity 
 % yeilds wavelength in channel widths measured down valley thru tortuosity
    waves(save_t).length=(diff(S_cum(abs(diff(temp))==1))/(B*tortuosity));
    if save_riv_vars == 1
        river(save_t).Xmig = dX_cl/dt; % migration distance
        river(save_t).Ymig = dY_cl/dt; % migration distance
    end
end
%% Compute new tortuosity and channel slope
% tortuosity calculated assuming straight valley
[~, tortuosity] = quicktor_TRH(X,Y);    
S_ch(t+1) = S_val/tortuosity;
%% Update flow variables
D_ch(t+1) = (Qo/W*sqrt(Cf_ch/g/S_ch(t+1)))^(2/3);
U_ch(t+1) = Qo/(D_ch(t+1)*W); 
% half-width:depth ratio for straight channel characterized by valley slope
%TRH (unused)   betaa(t+1) = B/D_ch(t+1); 
% Froude no^2 for straight channel characterized by valley slope
F2_ch(t+1) = (U_ch(t+1)/sqrt(g*D_ch(t+1)))^2;  
% clock_update = clock_update+toc; tic
%% Display progress.  Most of the code below here is written by TRH
if rem(t,disp_progress) == 0
    % TRH The river seems to extend off the downstream end of the valley 
    % over long simulation times, this truncates to original valley length
    % occasionally only at disp_progress interruptions
    new_end = find( X > x_max_orig , 1 ,'first');
    [new_end numel(X)]
    if new_end                       % not empty
        X(new_end:end)=[];           % remove the tail of the simulation
        Y(new_end:end)=[];
    end

    % Display progress:
    pct = t/sim_time*100;
    elapsed_algorithm_time = datenum(clock) - datenum(algorithm_start);
    disp([num2str(pct), ' % finished.',...
        sprintf(' Current Time:  %s \nest completion %s',...
        datestr(datenum(clock)) , datestr( datenum( algorithm_start ) +...
        elapsed_algorithm_time * 100 / pct ) ) ] )
    figure(1)
        plot(X,Y)
        title(sprintf('Final planform t is %d cfo is %f A is %d',...
            t,Cfo,alphaao))
        xlabel('Distance (m)'); ylabel('Distance (m)');
        axis equal
        drawnow
    figure(2)
        wavelength_straight = chan_len.*S_ch(1:length(bends))./...
            (bends*S_val*2*B);
        plot([1:length(bends)],wavelength_straight,[1:length(bends)],bends)
        title(sprintf(strcat('Number of bends, average average',...
            'wavelength %1.1f ± %0.2f of average'),...
            mean(wavelength_straight,'omitnan'),...
            std(wavelength_straight,'omitnan')))
        legend('Wavelength','Bends');
        xlabel('Model time steps'); 
        ylabel('Down valley AVERAGE meander wavelength, Bends');
        drawnow
    figure(3)
        plot(chan_len)
        title('channel length')
        xlabel('Model time steps'); ylabel('Length (m)');
        drawnow
    figure(4)
        py = freq_mig(1:end-1);     % end-1 because max value set by 
                                    % appending a number
        px = [1:length(py)]*mig_mag_max./length(py);
        bar(px,py)
        title('nodal migration rate using histogram not histc')
        xlabel(sprintf('Meters per %1.2f years',dt_yrs ) );
        ylabel('Frequency')
        drawnow
    figure(6)
        all_waves=vertcat(waves(:).length);
%         all_waves(1:10000)=[];  %suggest trimming off early measurements
        histogram(all_waves)
        title(sprintf('histogram of wavelengths avg %0.2f ± %0.3f',...
            mean(all_waves,'omitnan'),std(all_waves,'omitnan')))
        xlabel(sprintf('all waves Cfo %f, A %f, sim time %d years',...
            Cfo,alphaao,sim_time_yrs))
        ylabel('Frequency')
        drawnow
% DEBUG_migration_model_TRH_Ch2_line = 433
% keyboard

%% If all the timer variables and tic toc are uncommented, (above) then you 
% can display them with the code below including labels
%     a = [clock_curv clock_flow clock_angle clock_cutoff clock_smooth ...
%           clock_space clock_update];
%      timers = strcat('clock_curv clock_flow clock_angle clock_cutoff',...
%          ' clock_smooth clock_space clock_update total numel(riv)')
%     [a sum(a) numel(river)]

%     save('par_x_y.mat','X','Y','t','-v7.3'); %optional save output
%% code fragment for comparing different flowfield methods 'Error' is the 
% difference between methods
%     plot(ETM(ETM>0),EFM(EFM>0)); % comparing flowfield computations     
%     legend('Error truncated - mex','Error full - mex')
%     drawnow

end % if dislay progress
end % for t

%% finalize figures
figure(1)
    plot(X,Y)
    title(sprintf('Final planform %s',outfilename))
    xlabel('Down valley distance (m)')
    ylabel('Cross valley distance (m)')
    axis equal
%     print('-painters', '-dpng', '-r600', sprintf('%s planform.png',...
%         outfilename))
    %%
clock_time_end = clock;
end_time_is = num2str(clock_time_end)
duration = clock_time_end - clock_time_start

%% save run parameters
X=Xf; Y=Yf;         % taking penultimate coordinates which 
                    % match final ub calculation
% save(sprintf('%s_other.mat',outfilename),'scaleup','save_dt','dt_yrs',...
%     'sim_time_yrs','Cfo','alphaao','Eo','-v7.3')
% save(sprintf('%s_other_more.mat',outfilename),'scaleup','save_dt',...
%     'dt_yrs','sim_time_yrs','Cfo','alphaao','Eo','Cf_ch','S_val',...
%     'alphaa','So','Do','Uo','F2o','D_ch','S_ch','U_ch','F2_ch',...
%     'bends','all_waves','chan_len','X','Y','ub','-v7.3')
beep on % tell me you're finished
beep;pause(1);beep;
disp('You have reached the end of migration_model_TRH_Ch2.m')
disp('You have keyboard control of the variables.')
disp('You are encouraged to save variables and figures as desired')
disp('.svg and .fig are reccommended for editable figures.')
disp('Type dbcont or dbquit to end keyboard control')
keyboard % give user control dbquit or dbcontinue
end % function
