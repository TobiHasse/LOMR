function[Xreturn,Yreturn,numcuts,Xl,Yl,Xr,Yr,cutindices,ncut]=cutoff(Xs,...
    Ys,B,parallel)

% This function checks for cutoffs by seeing where the banks intersect
% themselves, and then performs the cut by removing the section of
% centerline. The point where the
% banks first intersect (looking downstream) is inserted into the
% centerline.
%
% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com
%
% INPUTS: Xs,Ys - centerline coordinate vectors
%         B - centerline half-width (scalar)
%
% OUTPUTS: Xreturn,Yreturn - new centerline coordinate vectors after
%                            cutoffs
%          numcuts - number of cutoffs performed in the time step
%          Xl,Yl and Xr,Yr - left and right bank coordinate vectors for
%                            post-cutoff centerline

%Xs=X;Ys=Y; 
% B=15; 

% Due to high curvatures (esp after a cutoff), the centerline can errantly 
% intersect itself. If this happens, cut off the loop. Eventually the 
% curvatures will smooth out.
cutindices = [];
ncut = [];

[~,~,icl,jcl] = intersections(Xs,Ys,1); % returns intersections of 
                                        % centerline with itself
while numel([icl; jcl]) > 0
    cutoff = floor(icl(1)); % assumes that only one i,j pair is returned 
                            % per intersection
    cuton = ceil(jcl(1));
    Xs(cutoff:cuton) = [];
    Ys(cutoff:cuton) = [];
    
    cutindices = [cutindices, cutoff, cuton];   % for statistics module 
                                                % and node tracking
    ncut = [ncut, cuton-cutoff+1];
    [~,~,icl,jcl] = intersections(Xs,Ys,1);
end

numsmooth = 10; % number of nodes before and after cutoff to use in 
                % smoothing spline
SGnfilt = 5;    % number of times to filter interpolated coordinates
SGorder = 3;    % order of the filter for interpolated 
                % coordinates -> lower order suppresses curvature more
add_pts = 15;   % number of points to add to interpolated section 
                % across cutoff

% Reshape inputs
Xs = Xs(:);
Ys = Ys(:);

% Return bank coordinates and intersection indices
if parallel
    [Xl,Yl,Xr,Yr,ir,jr,X0r,Y0r,il,jl,X0l,Y0l] = gen_banks_par(Xs,Ys,B);
else
    [Xl,Yl,Xr,Yr,ir,jr,X0r,Y0r,il,jl,X0l,Y0l] = gen_banks(Xs,Ys,B);
end

% Perform cutoffs
numcuts = 0;
while numel([ir;il]) > 0
    preXs = length(Xs);
    numcuts = numcuts + 1;
    % Cutoffs are performed moving downstream. This can make a difference
    % if there are multiple cutoffs in a time step.
    if numel(ir) == 0
              Xforcut = Xl;
              Yforcut = Yl;
              LR = 0;   % LR is just a marker that lets the code know 
                        % later if it's the right or left bank that was 
                        % cut; LR=0 ==> Left bank
        if numel(il)== 1
              LorRia = il(1); 
              LorRib = NaN; 
              LorRj = jl(1)-1;  % the -1 is to make up for the +1 added to 
                                % cuton later (exceeds length of stream if
                                %the intersection happens at the last node)
        else
              LorRia = il(1); 
              LorRib = il(2); 
              LorRj = jl(2);
        end

    elseif numel(il) == 0
              Xforcut = Xr;
              Yforcut = Yr;
              LR = 1;
       if numel(ir) == 1
              LorRia = ir(1);
              LorRib = NaN;
              LorRj = jr(1)-1;  % the -1 is to make up for the +1 added to 
                                % cuton later (exceeds length of stream if 
                                %the intersection happens at the last node)
       else
              LorRia = ir(1);
              LorRib = ir(2);
              LorRj = jr(2);
       end
    elseif min(ir) < min(il)
              LorRia = ir(1);
              LorRib = ir(2);
              LorRj = jr(2);
              Xforcut = Xr;
              Yforcut = Yr;
              LR = 1;
    else
              LorRia = il(1);
              LorRib = il(2);
              LorRj = jl(2);
              Xforcut = Xl;
              Yforcut = Yl;
              LR = 0;
    end
    
    if LorRib - LorRia == 1.5   % this lets the code know that it's dealing 
                                % with an end-of-the-river intersection, 
                                % pretty rare case but can happen. no 
                                % smoothing necessary
       cutoff = floor(LorRia);
       cutindices = [cutindices, cutoff, length(Xs)]; % for statistics 
                       % module and node tracking, the cuton is added later
       ncut = [ncut, length(Xs)-cutoff];
       Xs = Xs(1:cutoff);
       Ys = Ys(1:cutoff);
    else
        cutoff = floor(LorRia); %finds index at which banks intersect first
        cuton = ceil(LorRj)+1;  %finds index at which banks intersect last   

        cutindices = [cutindices, cutoff, cuton]; % for statistics module 
                            % and node tracking, the cuton is added later

        % Extract the centerline around the cutoff for interpolation and
        % smoothing

        numsmooth_beg = min(numsmooth, cutoff-1);   % in case the cutoff 
                                % occurs near the beginning of the river
        numsmooth_end = min(numsmooth, length(Xs)-cuton-1); % in case the 
                                % cutoff occurs near the end of the river

        Xsforsmooth = double([Xs(cutoff-numsmooth_beg:cutoff); ...
            Xs(cuton:cuton+numsmooth_end)]); % have to convert to double 
                                % to suppress a warning about mixed data 
                                % types from ode45 in interparc.
        Ysforsmooth = double([Ys(cutoff-numsmooth_beg:cutoff); ...
            Ys(cuton:cuton+numsmooth_end)]);
                
        % Interpolate across the cutoff
        [XY] = interparc(numsmooth*2+add_pts,Xsforsmooth,...
            Ysforsmooth,'spline');  % s_pts is a function of the 
                                    % streamlength, see further in code. 
                                    % The initial ratio of pts/length is 
                                    % maintained.
        Xsmooth = single(XY(:,1)); % convert back to a single 
        Ysmooth = single(XY(:,2));


        % Use Savitzky-Golay filter to smooth interpolated coordinates 
        % around cutoff

        window =  max(2.*round((SGorder+1+1)/2)-1,2.*...
            round((length(Xsforsmooth)/4+1)/2)-1); % window must be larger 
                            % than order of filter, hence the max(). 
                            % must also be odd, hence the rounding.
        tmp     =   [Xsmooth, Ysmooth];
        for iFilt=1:SGnfilt;
            tmp =   sgolayfilt(tmp,SGorder,window);
        end
        Xsmooth = tmp(:,1);
        Ysmooth = tmp(:,2);

        
        Xs = [Xs(1:cutoff-numsmooth-1);Xsmooth; Xs(cuton+numsmooth+1:end)];
        Ys = [Ys(1:cutoff-numsmooth-1);Ysmooth; Ys(cuton+numsmooth+1:end)];  
                
        % Get new banklines and cutoff indices after performing cutoff
        finalXs = length(Xs);
        ncut = [ncut, preXs - finalXs]; % for statistics module and node 
                                        % tracking
    end
    
    [Xl,Yl,Xr,Yr,ir,jr,X0r,Y0r,il,jl,X0l,Y0l] = gen_banks(Xs,Ys, B);
    
end

Xreturn = single(Xs);   % use single precision to save memory. in 
                        % spot-checks it made absolutely no detectable 
                        % difference but cuts file sizes (and therefore 
                        % memory use) in half
Yreturn = single(Ys);

end
