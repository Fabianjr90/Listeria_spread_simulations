function [beta] = fitAreaConvexHull(D_guess,k_guess,areaConvHullTime,areaConvHull)

lb = [0,0.1,1e-4];
ub = [];
options = optimset('Display','off');
fun = @(c,t)(c(1).*t.^2 - c(2).*t.*log(c(3).*t));
beta_guess = [pi*4*D_guess*k_guess,pi*4*D_guess,1];

beta = lsqcurvefit(fun,beta_guess,areaConvHullTime,areaConvHull,lb,ub,options);

end