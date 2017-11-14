function bactdif = chooseDcoeff(nbact,Dlow,Dhigh,conversionFactor)

% choose diffusion coefficient for nbact bacteria
randVals = rand(nbact,1);
bactdif = zeros(nbact,1);
bactdif(randVals <= conversionFactor) = Dlow;
bactdif(randVals > conversionFactor) = Dhigh;

end