function xyvals = levySample(D,alpha,n)
% sample from a polynomial-tailed distribution to generate a levy flight
% uses inverse tranform sampling

%random numbers from uniform distribution
uvals = rand(n,1);

%handle r<1 (if uvals < (alpha-1)/alpha))
rvals = uvals*alpha/(alpha-1); 

%handle r>1
ind2 = find(uvals>(alpha-1)/alpha);     
rvals(ind2) = (alpha*(1-uvals(ind2))).^(1/(1-alpha));

rvals = rvals*sqrt(D);

% pick directions randomly
thvals = rand(n,1)*2*pi;
xyvals(:,1) = rvals.*cos(thvals);
xyvals(:,2) = rvals.*sin(thvals); %sin^2 + cos^2 = 1
end