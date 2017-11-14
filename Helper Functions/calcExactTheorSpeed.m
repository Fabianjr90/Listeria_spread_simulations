function finalSpeed = calcExactTheorSpeed(D,k,time)

fun = @(c,t)( (-(4*D.*log(c(1)./(4*pi*D.*t))./t) - 4*D./t) ./ ...
    (2.*sqrt((4*D.*log(c(1)./(4*pi*D.*t)))./t + 4*D*k)) + ...
    sqrt((4*D.*log(c(1)./(4*pi*D.*t)))./t + 4*D*k));

FTG = 4.5;
myGuess = FTG*sqrt(D*time(end));
beta = 1/((1/(4*pi*D*time(end)))*exp(-(myGuess)^2/(4*D*time(end)) + k*time(end)));
finalSpeed = fun(beta,time(end));

end