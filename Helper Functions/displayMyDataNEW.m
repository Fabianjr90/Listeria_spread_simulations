function displayMyDataNEW(k_guess,D_guess,...
    observedSpeed,newObsSpeed,exactSpeed)

fprintf('\n')
display(['D_guess = ', num2str(D_guess)]);
display(['k_guess = ', num2str(k_guess)]);
fprintf('\n')

display(['c observed = ', num2str(observedSpeed)]);
display(['new c observed = ', num2str(newObsSpeed)]);
display(['exact speed = ',num2str(exactSpeed)]);
fprintf('\n')

cRatio = observedSpeed/exactSpeed;
newcRatio = newObsSpeed/exactSpeed;

display(['c ratio = ', num2str(cRatio)]);
display(['new c ratio = ', num2str(newcRatio)]);

end