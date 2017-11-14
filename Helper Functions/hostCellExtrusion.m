function [bactpos,nlive,nextruded] = hostCellExtrusion(bactpos,nlive,...
    hostCell,sizeHostCell)

% check which bacteria are within the host cell
bacteriaInHostCell = inpolygon(bactpos(1:nlive,1),...
    bactpos(1:nlive,2),hostCell(:,1),hostCell(:,2));
nextruded = sum(bacteriaInHostCell);

% plot host cell
hold on,plot(hostCell(:,1),hostCell(:,2),'g','LineWidth',4),hold off
pause(0.02)

% delete those bacterial entries
bactpos(bacteriaInHostCell,:) = [];
nlive = nlive-nextruded;

% calculate theta for all remaining bacteria
theta = atan2(bactpos(1:nlive,2),bactpos(1:nlive,1));

% move bacteria accordingly
bactpos(1:nlive,1) = bactpos(1:nlive,1)-sizeHostCell*cos(theta);
bactpos(1:nlive,2) = bactpos(1:nlive,2)-sizeHostCell*sin(theta);

end