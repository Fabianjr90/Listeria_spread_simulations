% Size of hop is sampled from normal distribution
% it measures radial speed

clc;clear;close all;
rng('shuffle');

%% Setting up simulation

nlive = 25;                           % initial number of bacteria
D = 1.0;                              % diffusion coefficient
delt = 0.01;                          % timestep
krep = 1.0;                           % replication rate
dTime = log(2)/krep;                  % doubling time
sigma = dTime/5;                      % width of distribution


maxnbact = 1e5;                       % maximum number of bacteria
maxstep = 5000;                       % maximum time
plotevery = 1;                        % how often you plot
plotwidth = 20;                       % for plotting 
com = [0,0];                          % center of mass

bactpos = zeros(maxnbact,2);          % positions of bacteria
bacthist = zeros(maxnbact,1);         % time experienced by bacteria
bactrep = zeros(maxnbact,1);          % doubling times of bacteria

% assign doubling times to initial nlive
bactrep(1:nlive) = normrnd(dTime,sigma,nlive,1);

% build fronts from the origin until finalFrontier
numFronts = 11;                      
finalFrontier = 500.0;
frontierPoints = repmat([0;1],numFronts,1)';
frontX = zeros(1,length(frontierPoints));
frontY = zeros(1,length(frontierPoints));
frontX((2:2:length(frontierPoints))) = finalFrontier.*cos(linspace(0,2*pi,numFronts));
frontY((2:2:length(frontierPoints))) = finalFrontier.*sin(linspace(0,2*pi,numFronts));

% Step 1: pre-defined vectors (for speed)
bacteria_number = zeros(maxstep,1);
avgR2 = zeros(maxstep,1);
time = zeros(maxstep,1);
time_convhull = zeros(maxstep,1);
area_convhull = zeros(maxstep,1);
front_distances = zeros(maxstep,1);

%% Simulation

% fullfig;
% profile on;
tic;

for step = 1:maxstep
        
    if (nlive>maxnbact)
        disp('Maximum bacteria reached')
        break
    end
    
    % calculate time
    time(step) = step*delt;
    bacthist(1:nlive) = bacthist(1:nlive)+delt;
    
    % Step 2: calculate size of hops (from Gaussian)
    delxy = randn(nlive,2)*sqrt(2*D*delt);
    bactpos(1:nlive,:) = bactpos(1:nlive,:)+delxy;
    
    % Step 3: calculate <r^2>
    dcom = bsxfun(@minus,bactpos(1:nlive,:),com);
    avgR2(step) = mean(sum(dcom.^2,2));
    
    % Step 4: calculate area of convex hull (plot every few steps)
    if mod(step,plotevery)==0
        [area_convhull(step),xp,yp] = ...
            simConvexHullNEW(bactpos(1:nlive,1),...
            bactpos(1:nlive,2),nlive,step,plotwidth);
        time_convhull(step) = delt*step;
        [xi,yi] = polyxpoly(frontX,frontY, xp, yp);
        pointInt = [xi(1:2:length(xi)),yi(1:2:length(yi))];
        front_distances(step) = mean(sqrt(sum(pointInt.^2,2)));  
    end
    
    % Step 5: calculate bacterial replication
    [bacthist,newbact,nrep] = bacteriaReplicationDT(...
        bactpos,bactrep,bacthist,nlive);
    bactpos(nlive+1:nlive+nrep,:) = newbact;
    bactrep(nlive+1:nlive+nrep,:) = normrnd(dTime,sigma,nrep,1);
    nlive = nlive+nrep;
    bacteria_number(step) = nlive;
    
end

toc;
% profile viewer;


%% Calculations

% calculating bacterial growth rate
g = fit(time(100:(step-1)),bacteria_number(100:(step-1)),'exp1');
figure,graphExpFit(time,bacteria_number,g,step);
k_guess = g.b;

% calculating slope of <r^2>
p = polyfit(time(100:(step-1)),avgR2(100:(step-1)),1);
figure, loglog(time,avgR2,'b.',time,time.*p(1) + p(2),'r-','linewidth',2);
graphLinearFit(time);
D_guess = p(1)/4;

% fitting area of convex hull
final = unique([time_convhull,area_convhull],'rows');
time_convhull = final(2:end,1);
area_convhull = final(2:end,2);
beta = fitAreaConvexHull(D_guess,k_guess,time_convhull,area_convhull);
area_fit = beta(1).*time_convhull.^2 -...
    beta(2).*time_convhull.*log(beta(3).*time_convhull);

% plotting convex hull data and fit
figure, plotConvexHull(final,time_convhull,area_fit);

% calculating speed of leading edge
r = (area_convhull.^0.5)/pi^0.5;
q = polyfit(time_convhull,r,1);
observedSpeed = q(1);
figure, plot(time_convhull,r,'b.',time_convhull(1:end),...
    time_convhull(1:end).*q(1)+q(2),'r-','linewidth',2);
graphLinearFitSpeed(time_convhull);

% calculating improved speed of leading edge
uniqueFront = unique(front_distances,'rows');
finalFront = uniqueFront(2:end);

s = polyfit(time_convhull,finalFront,1);
figure, plot(time_convhull,finalFront,'b.',time_convhull(1:end),...
    time_convhull(1:end).*s(1)+s(2),'r-','linewidth',2)
newObsSpeed = s(1);
graphLinearFitSpeed(time_convhull);
title('many fronts')

% calculate exact speed according to a reaction-diffusion equation
exactSpeed = calcExactTheorSpeed(D_guess,round(k_guess,1),time_convhull);

% displaying results
displayMyDataNEW(k_guess,D_guess,observedSpeed,newObsSpeed,exactSpeed);
