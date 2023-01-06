function [X, Y, interp_idcs] = enforce_spacing(X, Y,...
    dS_spacing_thresh, save_nodecount)
% PURPOSE: This function takes centerline coordinates (X,Y) and uses spline
%          interpolation to place nodes where the inter-node spacing 
%          surpasses the supplied threshold, dS_spacing_thresh. It works 
%          from upstream to downstream, interpolating one node at a time.
% AUTHOR:  Jon Schwenk, 2014. jonschwenk@gmail.com

% INPUTS: X, Y        - coordinates of channel centerline
%   dS_spacing_thresh - spacing beyond which a new node should be
%                       interpolated
%      save_nodecount - variable (either 1 or 0) that denotes whether 
%                       accounting should be done (for atom tracking)
% OUTPUTS: X,Y        - centerline coordinates with interpolated nodes
%         interp_idcs - indices of interpolated nodes: used for accounting
%                       algorithm

% keyboard

interp_idcs = []; % initialize

dS = sqrt(diff(X).^2+diff(Y).^2); % compute distance between nodes
insertnodes = find((dS > dS_spacing_thresh) == 1);  % find those that 
                                                    % surpass threshold

nloopcount = 0; % variables that helps detect problems 
nmax = numel(insertnodes);

preadd = 3;
postadd = 3;

while isempty(insertnodes) == 0
    % if statement to deal with cases where nodes are at the ends of the
    % centerline
    if insertnodes(1) < preadd+1
        preadd = insertnodes(1)-1;
        postadd = 3;
    elseif insertnodes(1) > length(X)-postadd-1
        preadd = 3;
        postadd = length(X)-insertnodes(1);
    end
    
    Xpiece = X(insertnodes(1)-preadd:insertnodes(1)+postadd);
    Ypiece = Y(insertnodes(1)-preadd:insertnodes(1)+postadd);
    
    tparam = 1:numel(Xpiece);
    insertptX = spline(tparam, Xpiece, ((preadd+1)+(preadd+2))/2); 
    insertptY = spline(tparam, Ypiece, ((preadd+1)+(preadd+2))/2);    
    
    X = [X(1:insertnodes(1)); insertptX; X(insertnodes+1:end)];
    Y = [Y(1:insertnodes(1)); insertptY; Y(insertnodes+1:end)];
    
    if save_nodecount
        interp_idcs = [interp_idcs; insertnodes(1)];
    end
    
    dS = sqrt(diff(X).^2+diff(Y).^2);
    insertnodes = find((dS > dS_spacing_thresh) == 1);
    
    nloopcount = nloopcount + 1;
    
    if nloopcount > nmax
       disp('enforce_spacing.m line 61 problem in node interpolation loop')
        keyboard
    end
end % while

end % function