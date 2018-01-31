% This script solves the steady state

% Intertemporal Euler equation
rss = 1/bett - (1-delt);
% Marginal product of capital
kss = (rss/alfa)^(1/(alfa-1));
% Production function
yss = ass * (muss * kss)^alfa * nss^(1-alfa);
% Marginal product of labor
wss = (1-alfa) * nss^(-alfa) * ass * (muss*kss)^alfa ;
% Capital law of motion
iss = delt * kss;
% Budget constraint
css = yss - iss;
% Marginal utility of consumption
lambdass = css^(-gamm);
% Intratemporal Euler
kappa = (wss * lambdass)/(nss^(sigg));

save parameterfile.mat;