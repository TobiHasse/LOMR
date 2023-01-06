function [C,cumlen,dS_btwn_nodes] = curvars_TRH(Xs, Ys, bc1, bcset, ...
    smoothnum)
% PURPOSE: This function is for computing the curvature between streamline 
%          centerline nodes, and also returns various other variables that 
%          are useful to other functions. 
%          Updates to the function 2015 include returning dS_btwn_nodes 
%          which is used by flowfield_TRH.  Previously dS (at nodes) was
%          returned but not used by anything.

% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com
%          Tobias Hasse, April 2015 tobiack@udel.edu

% INPUTS: Xs, Ys    - coordinates of channel centerline (or bankline)
%         bc1 bcset - curvature boundary conditions
%         smoothnum - curvature smoothing distance (in nodes)

% OUTPUTS: C        - curvature calculated via formulas of Johanesson and 
%                     Parker 1985 followed by Motta 2012
%         cumlen    - cumulative stream length in the streamwise coordinate 
%                     direction
%   dS_btwn_nodes   - used by flowfield

% % dCdS - dC/ds required by hydrodynamic model in calculation of a_p_5
% % del_s - stream length ds between nodes along stream centerline

% keyboard

% Xs=X;Ys=Y; % For testing only

% Reshape inputs
Xs = Xs(:);
Ys = Ys(:);

n_in = nargin;
if n_in == 2 % if only X,Y coordinates are given, set upstream C=0
    bcset = 2;
    smoothnum = 0;
end

%% Compute radius of circle fit through three consecutive points (other 
% methods coded at end of function)
x1 = Xs(1:end-2); y1 = Ys(1:end-2);
x2 = Xs(2:end-1); y2 = Ys(2:end-1);
x3 = Xs(3:end);   y3 = Ys(3:end);

At = 0.5 * ((x2-x1).*(y3-y1) - (y2-y1).*(x3-x1)); % signed areas of the 
                % triangle that has three consecutive points as vertices
sides_prod = (((x2-x1).^2+(y2-y1).^2).*((x3-x1).^2+(y3-y1).^2).*...
    ((x3-x2).^2+(y3-y2).^2)).^.5;
C = -4*At./sides_prod;

% Upstream boundary condition for curvature can be set different ways
if bcset == 1
    CS = bc1;   % 'periodic' boundary condition which imposes the 
                % curvature at the end of the river on the previous time 
                % step to the beginning of the river at the next
elseif bcset ==  2
    CS = 0; % C(1)=0. This causes the most-upstream portion of the river 
            % to straighten and elongate over time
elseif bcset == 3
    CS = mean(abs(C)) + std(abs(C))*randn; % random upstream curvature 
                                % (distribution has same mean and std as C)
else
    CS = 2*C(1)-C(2); % Do a linear interpolation. This seems to have the 
                      % same effect as as setting=0, but haven't tested it 
                      % over large times extensively.
end

% Downstream boundary condition for curvature is simply set to zero:
CE = 0;

C = vertcat(CS, C, CE);

%% Calculate stream length variables dS and cumlen

dS_btwn_nodes = sqrt(diff(Xs).^2+diff(Ys).^2); 
%TRH dS_at_nodes = .5*(dS_btwn_nodes(1:end-1)+dS_btwn_nodes(2:end));
% want dS at a node. this calculates it as the distance between a point 
% and the next point. therefore, average either side to get dS at node.
%TRH dS = [dS_at_nodes(1); dS_at_nodes; dS_at_nodes(end)]; 
% set dS boundaries equal to their nearest neighbors

cumlen = vertcat(0, cumsum(dS_btwn_nodes)); % this ensures that s=0 at 
                                            % the beginning of the curve

% %DEBUG
% d=diff(cumlen);
% e=(dS_btwn_nodes - d).^2;
% rmse = sqrt(sum(e)/numel(e));
% if rmse > 10^-10
%     rmse
% end

%TRH %% dcds by central differencing
% dC = diff(C);
% dCdS = dC./dS_btwn_nodes;
% dCdS = [dCdS; dCdS(end)];

%% Smooth Curvature
% Curvature smoothing as in Motta, eq'n 16, following Crosato 1990.
for u = 1:smoothnum
% for u = 1:5
    Cs = C(2:end-1);    % don't want boundary conditions to propogate
    C1 = Cs(1:end-2);
    C2 = Cs(2:end-1);
    C3 = Cs(3:end);
    Cs = vertcat(C(1), (2*Cs(1)+Cs(2))/3, (C1+2*C2+C3)/4, ...
        (2*Cs(end)+Cs(end-1))/3, C(end)); % reinsert boundary conditions
    C = Cs;
end

end

%% ----EXTRA CODE -- OTHER METHODS OF CURVATURE CALCULATION---%%
% Schwenk published an additionall 100 lines of code for this method, all
% commented out.  The interested coder is directed to Schwenk's 2015 paper
% and the supplementary material.  This code has been deleted by Tobias
% Hasse June, 2021

