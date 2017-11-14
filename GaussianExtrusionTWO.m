clc;clear;close all;
rng('shuffle');

%% Setting up simulation
% Listeria parameters

nlive = 1;                            % initial number of bacteria
orignlive = nlive;                    % const initial # bacteria
logGoal = 5;                          % orignlive x logGoal orders of mag
Dopt = 56.0;                          % optimal diffusion coefficient
Dother = 20.0;                        % "other" diffusion coefficient
delt = 0.01;                          % timestep
krep = 1.0;                           % replication rate
dTime = log(2)/krep;                  % doubling time
sigma = dTime/5;                      % width of distribution
cFactor = 0.95;                       % conversion factor

maxnbact = 6e5;                       % maximum number of bacteria
maxstep = 1e6;                        % maximum time
plotwidth = 100;                      % for plotting

% Step 1: pre-defined vectors (for speed)
cumExtrBact = zeros(1e6,1);      % number of extruded bacteria
bactpos = zeros(1e6,2);          % positions of bacteria
bacthist = zeros(1e6,1);         % time experienced by bacteria
bactrep = zeros(1e6,1);          % doubling times of bacteria
time = zeros(maxstep,1);         % total time in simulation
bactdif = zeros(1e6,1);          % diffusion coefficients

% assign doubling times to initial nlive
bactrep(1:nlive) = normrnd(dTime,sigma,nlive,1);
bactdif(1:nlive) = chooseDcoeff(nlive,Dopt,Dother,cFactor);

% host cell parameters
extrTime = 0.23;                                % how often hostCell extrudes
extrStep = extrTime/delt;                       % extrusion step size
extrSteps = [1,extrStep:extrStep:extrStep*1e5]; % all extrusion steps
extrCounter = 1;                                % count extrusion step
plotevery = extrStep/2;                         % how often you plot
sizeHostCell = 2.5;                             % radius of host cell
hostCell = createCircHostCell(sizeHostCell);    % create host cell

fullfig
tic;
for step = 2:maxstep
    
    % case 1: unfortunately, you have been wiped out    
    if (nlive==0)
        toc;
        disp('SORRY! All bacteria have died!')
        break
    
    % case 2: unfortunately, you have killed your host
    elseif (nlive>=maxnbact)
        plotExtrusion(bactpos(1:nlive,1),bactpos(1:nlive,2),...
            nlive,step,plotwidth,hostCell,...
            cumExtrBact(extrSteps(extrCounter)));
        toc;
        finalMessage = strcat(...
            'SORRY! The host died due to large bacterial load :(');
        disp(finalMessage)
        break
    
    % case 3: congrats! You extruded n logs your original number!
    elseif (cumExtrBact(step-1) > orignlive*10^logGoal)
        plotExtrusion(bactpos(1:nlive,1),bactpos(1:nlive,2),...
            nlive,step,plotwidth,hostCell,...
            cumExtrBact(extrSteps(extrCounter)));
        toc;
        finalMessage = strcat('CONGRATS! You have extruded ',{' '},...
            num2str(logGoal),' logs your initial number!');
        disp(finalMessage{1})
        break
    end
            
    % calculate time
    time(step) = step*delt;
    bacthist(1:nlive) = bacthist(1:nlive) + delt;
    cumExtrBact(step) = cumExtrBact(step-1);

    % Step 2: calculate size of hops (from heavy-tailed polynomial)
    delxy = randn(nlive,2).*sqrt(2*bactdif(1:nlive)*delt);
    bactpos(1:nlive,:) = bactpos(1:nlive,:)+delxy;
    
    % Step 3: calculate bacterial replication
    [bacthist,newbact,nrep] = bacteriaReplicationDT(...
        bactpos,bactrep,bacthist,nlive);
    bactpos(nlive+1:nlive+nrep,:) = newbact;
    bactrep(nlive+1:nlive+nrep,:) = normrnd(dTime,sigma,nrep,1);
    bactdif(nlive+1:nlive+nrep) = chooseDcoeff(nrep,Dopt,Dother,cFactor);
    nlive = nlive+nrep;
        
    % Step 4: plot and extrude bacteria
    if mod(step,plotevery)==0
        plotExtrusion(...
            bactpos(1:nlive,1),bactpos(1:nlive,2),...
            nlive,step,plotwidth,hostCell,...
            cumExtrBact(extrSteps(extrCounter)) );
        
        % extrude
        if (mod(time(step),extrTime)==0)
            [bactpos,nlive,extruBacteria] = ...
                hostCellExtrusion(bactpos,nlive,hostCell,sizeHostCell);
            
            % exit if you ran out of bacteria
            if (nlive==0)
                continue
            end
            
            % update how many bacteria have been extruded
            cumExtrBact(step) = cumExtrBact(step-1) + extruBacteria;
            
            % used to keep track of cumulative extruded bacteria
            extrCounter = extrCounter + 1;
            
            % update plot with bacteria that have been removed
            plotExtrusion(bactpos(1:nlive,1),bactpos(1:nlive,2),...
                nlive,step,plotwidth,hostCell,...
                cumExtrBact(extrSteps(extrCounter)));
        end
    end
end

% plot number of extruded bacteria vs time
if (nlive > 0)
    plotExtrudedBacteria(time(1:step-1),cumExtrBact(1:step-1),...
        extrSteps(1:extrCounter));
end
