function [Xl,Yl,Xr,Yr,ir,jr,X0r,Y0r,il,jl,X0l,Y0l] = gen_banks(Xsb,Ysb,B)
% This fucntion generates banks for a given centerline and constant width.
% INPUTS:   Xsb,Ysb - centerline coordinate vectors
%           B - channel half-width
% OUTPUTS:  Xl,Yl and Xr,Yr - left and right bank coordinate vectors
%           il,jl and ir,jr - left and right bank indices of where cutoffs
%                             occur
%           X0l,Y0l and X0r,Y0r - left and right bank intersection points 
%                                 corresponding to cutoffs found in i,j 
 
% There are two methods to generate banks. 

% METHOD 1: Uses sin(pi-angle)*B to
% calculate the new banks. This method underestimates the bank width
% everywhere, moreso at higher curvatures. However it's much faster and
% does a decent job so it's currently implemented.
% 

% METHOD 2: This code generates banks a distance B away from an input 
% centerline given by Xsb,Ysb. It does so by creating lines parallel to the
% vectors defining the centerline that are B distance away in a 
% perpendicular sense. The intersections of these parallel lines define the
% vetices of the banks. This method is currently pasted at the end of this 
% file.
%
% Both methods fail when curvature is high (e.g.immediately after a cutoff)
% causing the bank to erroneously intersect itself near the apex of the
% bend, sometimes multiple times. Extra coding was required to distinguish
% between an erroneous cutoff and a legitimate one.
%
% The code returns the indices of actual cutoffs and the left and right
% bank coordinates after erroneous bank intersections have been removed.

% Xsb=Xs;Ysb=Ys;  % for debugging/testing

% Reshape inputs
Xsb = Xsb(:);
Ysb = Ysb(:);

thresh = 7;
thresh2 = 10;
checkthresh = 5; % threshold for distinguishing between pairs of cutoff 
                 % vertices. if a whole section of bank is removed 
                 % (obviously) incorrectly, increase this threshold
min_ox_length = 20;
% n_gen_line = 10^5;    % just needs to be large enough so that the line 
                        % segments overlap

% angles computed under assumption of constant valley direction found via 
% linear regression!!!
[angles] = cl_angles(Xsb,Ysb);  

Xl = Xsb-sin(pi-angles)*B;
Yl = Ysb-cos(pi-angles)*B; 
Xr = Xsb+sin(pi-angles)*B;
Yr = Ysb+cos(pi-angles)*B; 

% Loop to generate left and right banks
for lr = 1:2
    if lr == 1; % Determines left or right bank
        Xbank = Xr;
        Ybank = Yr;
        B = B; % right bank
    else
        Xbank = Xl;
        Ybank = Yl;
        B = -B; % left bank
    end    
    
    % Now that the bank has been generated, need to take care of erroneous
    % cases.
    % returns intersections of left bank with itself
    [X0fix,Y0fix,ifix,jfix] = intersections(Xbank,Ybank,1); 
    
    % This if statment checks for different combinations of legitmate
    % cutoffs mixed in with erroneous bank intersections
    if isempty(ifix) == 1
    elseif length(ifix) == 1 && (length(Xbank) - floor(jfix)) < 3 % case 
        % where the end of the river is intersecting (only one intersection
        % point in this case, hence the isodd requirement)
        ifix = [ifix ifix+1.5]; % adding 1.5 as a marker to the cutoff.m 
                                % routine that this is an endpoint
        jfix = [jfix-1.5 jfix]; % adding 1.5 as a marker to the cutoff.m 
                                % routine that this is an endpoint
        X0fix = [X0fix X0fix];
        Y0fix = [Y0fix Y0fix];
        
    elseif length(ifix) == 1 % case where there is only one bank 
        % intersection (must be a bank error in this case UNLESS the 
        % cutoff points intersect EXACTLY at the same node, extremely 
        % low likelihood)
        cutoff = floor(ifix);
        cuton = ceil(jfix);
        Xbef = Xbank(1:cutoff);
        Xaft = Xbank(cuton:end);
        Xbank = [Xbef; X0fix; Xaft];
        Ybef = Ybank(1:cutoff);
        Yaft = Ybank(cuton:end);
        Ybank = [Ybef; Y0fix; Yaft];
        ifix = [];
        jfix = [];
              
    else    % case where valid cutoffs are mixed in with SINGLE CROSS bank 
        % error intersection
        keep = find(jfix-ifix > min_ox_length);   % if the cutoff section 
        % is fewer than min_ox_length number of vertices, it's erroneous
        ijXY = [ifix jfix X0fix Y0fix];   % puts all the intersections in 
                                          % matrix form for less coding
        remmat = ijXY;  % makes a copy of the intersection matrix
        remmat(keep,:) = [];    % leaves only the intersections that 
                                % should be removed
        keepmat = ijXY(keep,:); % a matrix of the intersections to keep; 
                                % each row corresponds to an intersection. 
                                % should be even number of rows!
        
        ifix = keepmat(:,1);    % leaves all the coupled intersection 
                                % indices that are actual cutoffs
        jfix = keepmat(:,2);
        X0fix = keepmat(:,3);
        Y0fix = keepmat(:,4);
        
        icheck = floor(remmat(:,1));
        jcheck = ceil(remmat(:,2));
        % These loops are required to handle cases where the returned
        % intersections are contained within one another; that is, an
        % intersection may be at i=20 and j=30, but there's another one at
        % i=22 and j=27 and another one at i = 19 and j = 29, for example.
        for xx = 1:length(icheck)
            for intno = 1:size(remmat,1)
                if icheck(xx) > floor(remmat(intno,1)) && icheck(xx) < ...
                        ceil(remmat(intno,2))
                end
                if jcheck(xx) > floor(remmat(intno,1)) && jcheck(xx) < ...
                        ceil(remmat(intno,2))
                end
            end
        end
        
        for g = 1:length(icheck)-1
            if icheck(g+1) == 0 && jcheck(g+1) == 0
            elseif icheck(g+1) == 0
                if jcheck(g+1) > jcheck(g)&& abs(jcheck(g+1)-jcheck(g))<...
                        checkthresh
                    jcheck(g) = jcheck(g+1);
                end
                jcheck(g+1) = 0;
            elseif jcheck(g+1) == 0
                if icheck(g+1) < icheck(g)&& abs(icheck(g+1)-icheck(g))<...
                        checkthresh
                    icheck(g) = icheck(g+1);
                end
                icheck(g+1) = 0;
            end
        end
        ijcheck = [icheck jcheck];
        ijcheck(icheck==0 | jcheck==0,:)=[];
        icheck = ijcheck(:,1);
        jcheck = ijcheck(:,2);
        if length(icheck)~=length(jcheck)
            disp('error in check vectors, gen banks')
        end
       
        removeindices = [];
        for r = 1:length(icheck)
            removeindices = [removeindices icheck(r):jcheck(r)];
        end
        removeindices = unique(removeindices);
        
        % Now that the indices to remove have been found, remove them.
        Xbank(removeindices) = [];
        Ybank(removeindices) = [];
    end
    
    % Save the coordinates and intersections
    if lr == 1
        Xr = Xbank;
        Yr = Ybank;
        ir = ifix;
        jr = jfix;
        X0r = X0fix;
        Y0r = Y0fix;
    else
        Xl = Xbank;
        Yl = Ybank;
        il = ifix;
        jl = jfix;
        X0l = X0fix;
        Y0l = Y0fix;
    end
    
    clear Xbank Ybank res ifix jfix X0fix Y0fix cutout n_vertices_btwn

end % l/r end

end % function end

%% Method 2  
% Schwenk published an additionall 200 lines of code for this method, all
% commented out.  The interested coder is directed to Schwenk's 2015 paper
% and the supplementary material.  This code has been deleted by Tobias
% Hasse June, 2021
