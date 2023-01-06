function save_params_deposition()
% Purpose:  organize many of the adjustable inputs for sediment deposition
%           and converting from a node based channel centerline to a pixel
%           based channel and floodplain
%           Create this file in commented code for better communication and
%           understanding
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     June, 2021

load params_storage.mat
load params_meander.mat
% see equation 1 from Tobias Hasse's dissertation, the floodplain sediment
% thickness ranges from 0 to 1
pt_bar_elev = .35; % min depth of deposit after channel abandons pixel
lambda = .5;       % diffusive exponent
% Alan Howard (1992) used these values in his sedimentation equation
% Mu = .4;         % diffussive sed rate (coarse sediment) per 100 years
% Nu = .125;       % 'pelagic' sedimentation all over (fine sed) / 100 yr
Mu = .4   * time_step_years/100;    % diffussive sed rate per time step
Nu = .125 * time_step_years/100;    % 'pelagic' sed rate per time step
% NOTE: 'mu' is a MATLAB function, so Mu and Nu are used as variables

px_size_thresh = 100; % see rm_islands.m for discussion of this threshold
% depends on migration rate, meander wavelength, pixel size, dt

pixel_sz = 2 * B / pix_per_chan;  % size of pixels in meters

% interpolated nodes per pixel when converting centerline nodes to raster
% with too few nodes, the single pixel channel centerline has gaps, with
% too many nodes it takes extra computation.
interp_nds_pix = 3;     % 3 is a good number                 

save('params_deposition.mat','pt_bar_elev','lambda','Mu','Nu',...
    'px_size_thresh','pixel_sz','interp_nds_pix','-v7.3')

end