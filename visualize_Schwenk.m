% Purpose:  visualize centerline node simulation results
%           This will play an animation of the entire simulation moving
%           forward through time, then, if "nodecount" is available it will
%           step through each cutoff meander bend and play an animation
%           moving backward through time from the time of cutoff to the
%           'beginning' of the channel segment which becomes the cutoff
%           bend
% Authors:   Schwenk, 2014, lightly edited by Hasse 2015, 2021

%% Visualize Schwenk model results (centerline migration)
% Find plot limits
figure(1)
n = numel(river);
pstep = 5; % if save_dt = 20, pstep = 5 is 10 years (dt = .1 years)
Xmin = min(vertcat(river(1:pstep:n).Xcl)); Ymin = min(vertcat(river(1:pstep:n).Ycl));
Xmax = max(vertcat(river(1:pstep:n).Xcl)); Ymax = max(vertcat(river(1:pstep:n).Ycl));
windowadd = 500;
%% Plot
for p = 1:pstep:n
    clf
    plot(river(p).Xcl,river(p).Ycl,'marker','o','markersize',1,'linestyle','none')
    hold on
    xlim([Xmin-windowadd Xmax+windowadd]); 
    ylim([Ymin-windowadd Ymax+windowadd]);
    axis equal
    if p == 1 
        disp('press any key to play movie')
        title('press any key to play movie')
        waitforbuttonpress
    end
    title(['t = ',num2str(p)]);
    pause(.05);
end

if isequal(nodecount,0)
    fprintf(['nodecount is empty.  Try re-running the simulation with'...
        ' save_nodecount = 1. \n Change this in migration_model_',...
        'TRH_Ch2 at line 66\n'])
    return
end
%% Extract individual atoms through time (ss)
atoms = atom_tracking_TRH(river, nodecount, B);

%% Visualize individual bend evolutions
backwards = 1; % 0 ends in cutoff, 1 starts with cutoff and evolves backwards
group = 1:numel(atoms);
nskip = 10;
for j =1:numel(group)
    mID = group(j);
    plotend = numel(atoms(mID).t);
    tend = atoms(mID).t(1);
    t1 = atoms(mID).t(end);
    idx1end = atoms(mID).idx1(1);
    idx2end = atoms(mID).idx2(1);
    idx11 = atoms(mID).idx1(end);
    idx21 = atoms(mID).idx2(end);
    
    windowadd = 500;
    xlimsearch = [river(t1).Xcl(idx11:idx21); river(tend).Xcl(idx1end:idx2end)];
    ylimsearch = [river(t1).Ycl(idx11:idx21); river(tend).Ycl(idx1end:idx2end)];
        xlimit = [min(xlimsearch)-windowadd max(xlimsearch)+windowadd]; 
        ylimit = [min(ylimsearch)-windowadd max(ylimsearch)+windowadd]; 
                 
        if backwards == 1
            fb = 1:nskip:plotend;
        else
            fb = plotend:-nskip:1;
        end
        for q = 1:floor(length(fb)/1)
            p = fb(q);
            X_mdr = river(atoms(mID).t(p)).Xcl(atoms(mID).idx1(p):atoms(mID).idx2(p));
            Y_mdr = river(atoms(mID).t(p)).Ycl(atoms(mID).idx1(p):atoms(mID).idx2(p));

            X = river(atoms(mID).t(p)).Xcl;
            Y = river(atoms(mID).t(p)).Ycl;

            plot(X,Y,'b.'); hold on
            plot(X_mdr,Y_mdr,'r.')
            axis equal
            xlim(xlimit); ylim(ylimit)
            
            str=sprintf('Bend number %d; t = %d', mID, atoms(mID).t(p));
            title(str)

            if q == 1 
                waitforbuttonpress
            end
            pause(0.007)     
            clf
        end
    if j ~= numel(group)
        disp('Press any key or click mouse for next bend...')
        waitforbuttonpress
    end
   
end
