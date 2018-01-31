// ----------------
// Utlilization cost model 
// Author: Anh Trinh
// Date: November, 2017
// ----------------

// ----------------
// Endogenous variables
// ----------------
var
    ct, yt, nt, at, it, ktp1, wt, rt, mut, lambdat;

// ----------------
// Exogeneous variables
// ----------------
varexo
    etaat;

// ----------------
// Parameters
// ----------------
parameters
    alfa, bett, delt, rho, varsigmaa, gamm, curvv, sigg, ass, rss, kappa ;

// ----------------
// Calibrated parameters
// ----------------
// parameterfile.mat is a MAT file where we had stored the parameter values
// and the steady-state values.
load parameterfile;
set_param_value('alfa',alfa);
set_param_value('bett',bett);
set_param_value('delt',delt);
set_param_value('rho',rho);
set_param_value('varsigmaa',varsigmaa);
set_param_value('gamm',gamm);
set_param_value('curvv',curvv);
set_param_value('sigg',sigg);
set_param_value('ass',ass);
set_param_value('rss',rss);
set_param_value('kappa',kappa);

// ----------------
// Model step
// ----------------
model;
// 1. Productivity process
    log(exp(at)) = (1-rho) * log(ass) + rho * log(exp(at(-1))) + varsigmaa * etaat;
// 2. Capital law of motion
    exp(ktp1) - (1-delt) * exp(ktp1(-1)) = exp(it);
// 3. Marginal utility
    exp(ct)^(-gamm) = exp(lambdat) ;
// 4. Intratemporal Euler equation
    exp(lambdat) * exp(wt) = exp(nt)^sigg * kappa;
// 5. Intertemporal Euler equation
    exp(ct)^(-gamm) = bett * exp(ct(+1))^(-gamm) * (exp(rt(+1)) * exp(mut(+1)) - rss * ((exp( curvv * (exp(mut(+1))-1) )-1)/curvv) + 1 - delt ) ; 
// 6. Production function
    exp(yt) = exp(at) * (exp(mut) * exp(ktp1(-1)))^alfa * exp(nt)^(1-alfa);
// 7. Marginal product of capital
    exp(rt) = alfa * exp(yt)/ (exp(mut) * exp(ktp1(-1))) ;
// 8. Marginal product of labor
    exp(wt) = (1-alfa) * exp(yt) / exp(nt);
// 9. Resource constraint
    exp(yt) = exp(ct) + exp(it) + rss * ( (exp( curvv * (exp(mut)-1) )-1)/curvv );
//10. FOC of capacity utilization
    exp(rt) = rss * exp((exp(mut) - 1)* curvv);
end;

// ----------------
// Steady state step
// ----------------
initval;
    ct = log(css);
    yt = log(yss);
    nt = log(nss);
    it = log(iss);
    ktp1 = log(kss);
    at = log(ass);
    wt = log(wss);
    rt = log(rss);
    lambdat = log(lambdass);
    mut = log(muss);
end;

steady;
check;

// ----------------
// Shock step
// ----------------
shocks;
    var etaat; stderr 1;
end;

// ----------------
// Stochastic simulation step
// ----------------
stoch_simul(order=1,relative_irf,irf=40,pruning);
