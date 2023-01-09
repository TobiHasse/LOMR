function [ub] = flowfield_TRH( B, C, Do, s, dS, alphaa, F2, Cf_ch, ...
    truncate_conv, node_thresh)
% 11/19/2014 - Did speed tests on for computing integral two additional
% ways: 1) over-allocating the integrands vector, padding the extras with
% zeros so that Matlab doesn't need to keep creating a new one each
% iteration. 2) Computing the convolution in a for loop and truncating when
% the additional contribution falls below some threshold. Both methods were
% longer. The for loop was exceptionally longer.

% Updated by Tobias TRH March/April 2015 tobiack@udel.edu
% took out test if sidx == 1: do operation before loop, loop 
% from 2:numel(s) vectorized us_idx using function near_elem_TRH, major 
% speed improvement ±25. took dS computation out of loop, pass it in from
% curvar.m access by indexing into vector using 'end' method of vector 
% indexing flowfield_TRH more than 2x faster than flowfield_Shwenk 

% compute parameters and nondimensionalize variables
s = s/B;               % downstream distance, normalized by half width
dS= dS/B;              % normalize dS from curvars_TRH.m by half width
Ro = min(abs(1./C));   % minimum radius of curvature along reach
C = C*Ro; 
beta = B/Do;           % (half-width) aspect ratio  
%TRH: using current channel depth rather than straightened depth tempers 
% the change in lam_o requiring larger changes in Cfo to observe an effect
lam_o = -2*beta*Cf_ch; % characteristic expon describing upstream influence
nu_o = B/Ro;

% ************************************************************************
A = alphaa + 1;%accounts for secondary flow, see Sun 1996 among many others 
% TRH the +1 apparently comes from JP '85, equation 14 & 15 but the 
% reference is in ERROR!!!!!!!!!!!!!! and it should be A = alphaa - 1
% ************************************************************************

N = numel(s);
% dS is already created in curvar to compute cumlen (s), just pass it in
% dS = diff(s);          %TRH compute dS out here rather than in loop

% TRH optimizing cutoff search, how short can a cutoff be? 
% see fig 7 of Schwenk
% Lo_in_Bs = 1/abs(lam_o);

% upstream boundary condition
ub(1) = 0;
if truncate_conv 
%     % Truncate convolution integral at thresh distance upstream (or at 
%     % source node)
%     conv_int_thresh = 100; % multiples of B: threshold distance to 
%     % compute convolution integral (saves computation time) 
%     % conv_int_thresh is made obsolete by node_thresh
%     conv_int_thresh = node_thresh;

%     % Vectorize computation of the upstream index in a function call
%     us_idx = near_elem_TRH(s-conv_int_thresh,s); % 40% faster 

%     % Ideally flow_conv_int_trunc could recieve us_idx from above, 
%     % however I cant figure out how to pass a vector of indexes to 
%     % a *.c (mex) function.
%     % The node_thresh approximates the maximum length of the truncated
%     % convelution integral.  This mex is 17% faster than the loop below
    try
        int_term = flow_conv_int_trunc ( N, node_thresh, lam_o, C, s, dS );
    catch
        fprintf(['Error running c code in function: \n"%s" \n',...
            'see comments in file \n'],...
            fullfile(fileparts(mfilename('fullpath')),mfilename))
        keyboard
    end
%   % If you have not tried to compile the c function try typing:
%   % mex flow_conv_int_trunc.c  (on the command line)
%   % if that fails UNCOMMENT the next 8 lines (and us_idx Line 49 above)
%   % to compute the int_term in this function not in the *.c function call
%   int_term = NaN(N,1);
%   int_term(1)=0;
%   % Depending on channel length, truncating the convelution interval
%   % in this loop might be slower than not truncating it (Line 90 below)
%     for ds_idx = 2:N 
%             up_id = us_idx(ds_idx);
%             integrands = C(up_id:ds_idx).*exp(lam_o*(s(ds_idx)-...
%                 s(up_id:ds_idx)));
%             int_term(ds_idx) = dS(up_id:ds_idx-1)' * ...
%                 (integrands(1:end-1) + integrands(2:end))/2;
%     end
%% for comparing methods of computing the int_term
%     % change one of the int_term variables above to int_t
%     if max((int_term - int_t).^2) > 10^-15  % test for accuracy
%         sprintf('new method is not computing int term right')
%         max((int_term - int_t).^2)
%         keyboard
%     end


else % compute convelution integral without truncation
    % compiled function flow_conv_int (*.c mex) is slowest method (OS & 
    % compiler specific) suggest truncating (above) or running loop below 
%     int_term =   flow_conv_int ( N, lam_o, C, s, dS);

%     % compute convolution integral in loop
    int_term = NaN(N,1);
    int_term(1)=0;
    for sidx = 2:N
            integrands = C(1:sidx).*exp(lam_o*(s(sidx)-s(1:sidx)));
            int_term(sidx) = dS(1:sidx-1)' * (integrands(1:end-1) +...
                integrands(2:end))/2;
    end
%% for comparing methods of computing the int_term
%     % change one of the int_term variables above to int_t
%     if max((int_term - int_t).^2) > 10^-9  % test for accuracy
%         sprintf('new method is not computing int term right')
%         max((int_term - int_t).^2)
%         keyboard
%     end
    
%% for comparing methods of computing the int_term
% could also taking the abs() of each term before comparing?
% etm = sum( ( int_term      - int_term_mex ).^2 );
% efm = sum( ( int_term_full - int_term_mex ).^2 );
% eft = sum( ( int_term_full - int_term     ).^2 );
% 
%     if (etm>10^-6 & efm>10^-6)
%         if (eft < efm)
%             sprintf(strcat('mexfile flow_conv_int.c is not computing',...
%                 ' int term right'))
%             eft - efm
%         end
%     end

end
% equation 2 in Schwenk 2015
ub = (ub(1)+nu_o*(C(1))).*exp(lam_o*s) + nu_o.* ...
    (-C - lam_o/2.*(F2 + A).*int_term); 

end % function

%% EXTRA CODE
% cut from below computation of us_idx:
% m and preallocation of integrands may be useful if subsequent loop
% is in a compiled subfunction
%     m = [1:N]-us_idx+1; % TRH the max of this would give the longest 
%                         % integrands vector
%     integrands_length = max(m);   % TRH preallocation perhaps for 
%                                   % call to compiled code
%     integrands = zeros(integrands_length,1);
