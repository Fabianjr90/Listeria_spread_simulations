function [bacthist,newbact,nrep] = bacteriaReplicationDT(...
    bactpos,bactrep,bacthist,nlive)

repind = find(bactrep(1:nlive)<=bacthist(1:nlive));
nrep = length(repind);
newbact = bactpos(repind,:);
bacthist(repind) = 0;              %reset history of replicated bacteria

end