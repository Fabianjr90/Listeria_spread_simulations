% bacteria will come into contact with boundaries
% there is gamma probability of them crossing boundary
% there is gamma-1 probability of staying at the boundary

clc;clear;close all;
rng('shuffle');

%% Setting up simulation

nlive = 5;                           % initial number of bacteria
D = 0.1;                              % diffusion coefficient
delt = 0.01;                          % timestep
krep = 1.0;                           % replication rate
dTime = log(2)/krep;                  % doubling time
sigma = dTime/5;                      % width of distribution
gamma = 0.0;                          % probability of crossing boundary
beta = 0.3;

maxnbact = 1e5;                       % maximum number of bacteria
maxstep = 5000;                       % maximum time
plotevery = 1;                        % how often you plot
plotwidth = 1.6;                       % for plotting 

bactangle = zeros(maxnbact,1);        % last angle
bactpos = zeros(maxnbact,2);          % positions of bacteria
bacthist = zeros(maxnbact,1);         % time experienced by bacteria
bactrep = zeros(maxnbact,1);          % doubling times of bacteria
bactrep(1:nlive) = ...                % assign doubling times to nlive
    normrnd(dTime,sigma,nlive,1);

% pre-defined vectors (for speed)
bacteria_number = zeros(maxstep,1);
time = zeros(maxstep,1);

% assign initial random angles (uniform)
bactangle(1:nlive) = rand(nlive,1)*2*pi;

% host cell boundaries and coordinates
% N = 12.5;
N = 12.5;
XY_h = [ -N*ones(2*N+1,1), (-N:1:N)', N*ones(N*2+1,1),(-N:1:N)' ];
XY_v = [ (-N:1:N)', -N*ones(2*N+1,1), (-N:1:N)', N*ones(2*N+1,1)];
boundaries = [XY_h;XY_v];
boundLength = length(boundaries);
boundCoord = zeros((N+1)*4-1,2);
boundCoord(2:length(XY_h)+1,:) = [zeros(length(XY_v),1),XY_v(:,1)];
boundCoord(length(XY_v)+2:end,:) = [XY_h(:,2),zeros(length(XY_h),1)];

% alternative boundary plotting
[X,Y] = meshgrid(-N:1:N,-N:1:N);

%% Simulation

% fullfig;
% profile on;
tic;

for step = 1:maxstep
        
    if (nlive>maxnbact)
        disp('Maximum bacteria reached')
        break
    end
    
    % step 1: calculate time
    time(step) = step*delt;
    bacthist(1:nlive) = bacthist(1:nlive)+delt;
    
    % Step 2: calculate size of hops (from Gaussian)
    thvals = normrnd(bactangle(1:nlive),beta);
    bactangle(1:nlive) = thvals;
    delxy = [cos(thvals)*sqrt(2*D*delt),sin(thvals)*sqrt(2*D*delt)];
    temp = bactpos(1:nlive,:)+delxy;
        
    % Step 3: check for bouncers (cross if mygamma is less than gamma)
    out = lineSegmentIntersect([bactpos(1:nlive,:),temp],boundaries);
    allInts = out.intAdjacencyMatrix;
    tocross = sum(allInts,2)>0;
    mygammas = rand(nlive,1).*tocross;    
    bouncer = mygammas > gamma;
    
    % Step 4: calculate distances from origin to intersections
    intX = out.intMatrixX(bouncer,:);
    intY = out.intMatrixY(bouncer,:);
    distances = sqrt( (intX-bactpos(bouncer,1)).^2 + ...
        (intY-bactpos(bouncer,2)).^2) + 100*D.*~out.intAdjacencyMatrix(bouncer,:);
    upperLimit = max(distances(:));
    [minDist, whichbound] = min(distances,[],2);
    trueBound = whichbound;
    whichbound = whichbound+1;
    
    % If it was supposed to cross two membranes, and succeeded in the
    % first, then bounce off the second.
    dCrossers = sum(allInts>0,2)>1;
    dCrossersFS = logical(dCrossers.*~bouncer);
    
    if sum(dCrossersFS)>0
        xFS = out.intMatrixX(dCrossersFS,:);
        yFS = out.intMatrixY(dCrossersFS,:);
        distancesFS = sqrt( (xFS-bactpos(dCrossersFS,1)).^2 + ...
            (yFS-bactpos(dCrossersFS,2)).^2).*out.intAdjacencyMatrix(dCrossersFS,:);
        [~, whichboundFS] = max(distancesFS,[],2);
        whichboundFS = whichboundFS + 1;
        nFS = bactpos(dCrossersFS,:).*(boundCoord(whichboundFS,:)~=0) - ...
            boundCoord(whichboundFS,:);
        bactpos(dCrossersFS,:) = bactpos(dCrossersFS,:) - 2*nFS;
    end
    
    % Step 5: move bacteria appropriately
    bactpos(1:nlive,:) = bactpos(1:nlive,:) + delxy;
    n = temp(bouncer,:).*(boundCoord(whichbound,:) ~= 0)-boundCoord(whichbound,:);
    bactpos(bouncer,:) = bactpos(bouncer,:) - 2*n;
        
    % Step 6: check for second bounce
    [nrows,~] = size(intX);
    linearInds = nrows*(trueBound-1) + (1:nrows)';
    newOrigins = [intX(linearInds),intY(linearInds)];
    newOut = lineSegmentIntersect(...
        [newOrigins,bactpos(bouncer,:)],boundaries);
    isDoubleBouncer = sum(newOut.intAdjacencyMatrix,2) > 1;
    
    % Step 7: bounce the double bouncers
    doubleBouncer = zeros(length(bouncer),1);
    doubleBouncer(bouncer) = isDoubleBouncer;
    doubleBouncer = logical(doubleBouncer);
    doubleX = newOut.intMatrixX;
    doubleY = newOut.intMatrixY;
    dDistances = sqrt( (doubleX - newOrigins(:,1)).^2 + ...
        (doubleY - newOrigins(:,2)).^2).*...
        newOut.intAdjacencyMatrix;
    [maxDist, n_whichbound] = max(dDistances,[],2);
    truenbound = n_whichbound;
    n_whichbound=n_whichbound.*isDoubleBouncer+1;
    n_n = bactpos(bouncer,:).*(boundCoord(n_whichbound,:) ~= 0)...
        -boundCoord(n_whichbound,:);
    bactpos(bouncer,:) = bactpos(bouncer,:) - 2.*n_n;
    
    % Step 8: recalculate angle
    newthvals = atan2( (bactpos(bouncer,2) - newOrigins(:,2) ),...
        (bactpos(bouncer,1) - newOrigins(:,1) ));
    bactangle(bouncer) = newthvals;

    % Step 9: calculate bacterial replication
    [bacthist,newbact,nrep] = bacteriaReplicationDT(...
        bactpos,bactrep,bacthist,nlive);
    bactpos(nlive+1:nlive+nrep,:) = newbact;
    bactrep(nlive+1:nlive+nrep,:) = normrnd(dTime,sigma,nrep,1);
    bactangle(nlive+1:nlive+nrep) = rand(nrep,1)*2*pi;
    nlive = nlive+nrep;
    bacteria_number(step) = nlive;
    
    % Step 10: optional plotting
    if mod(step,plotevery)==0
        plot(X,Y,'Color',[0.8 0.8 0.8],'LineWidth',2)
        hold on
        plot(Y,X,'Color',[0.8 0.8 0.8],'LineWidth',2)
        scatter(bactpos(1:nlive,1),bactpos(1:nlive,2))
%         gscatter(bactpos(1:nlive,1),bactpos(1:nlive,2),1:nlive,'','',30)
%         scatter(bactpos(1:nlive,1),bactpos(1:nlive,2),300,'r.')
        legend('off')
        set(gca,'FontSize',20)
        title(sprintf('Step %d, nlive: %d',step, nlive))
        hold off
        pbaspect([1 1 1]);
        xlim([-plotwidth,plotwidth])
        ylim([-plotwidth,plotwidth])
        drawnow
    end
    
end
toc;
% profile viewer;

x = bactpos(:,1);
y = bactpos(:,2);

[k,~] = boundary(x,y,1);
data = [y(k),x(k)]';
[z, r, L2norm, residuals] = fitcircle(data);
MSE = 1/length(residuals) * sum(residuals.^2);

scatter(y,x,'.k')
hold on;
plot(y(k),x(k),'r','LineWidth',3);
t = linspace(0, 2*pi, 100);
plot(z(1)  + r  * cos(t), z(2)  + r * sin(t), 'b','LineWidth',3)
axis equal

% figure,plot(X,Y,'Color', [0.8 0.8 0.8])
% hold on
% plot(Y,X,'Color', [0.8 0.8 0.8])
% scatter(bactpos(1:nlive,1),bactpos(1:nlive,2),100,'r.')
% xlim([-10 10])
% ylim([-10 10])
% pbaspect([1 1 1]);
