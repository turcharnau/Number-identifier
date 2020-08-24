clear;
%
% Parameters for dataset generation
%
num_target =[8];
tr_freq    = .5;   %tots els números amb la mateixa frequència      
tr_p       = 250;       
te_q       = 250;       
tr_seed    = 369354;    
te_seed    = 169306;    
%


% Parameters for optimization
%
la = 10.0;                                                     % L2 regularization.
epsG = 1.0e-6; kmax = 10000;                                   % Stopping criterium.
ils=1; ialmax = 2; kmaxBLS=30; epsal=1.0e-3;c1=0.01; c2=0.45;  % Linesearch.
isd = 7; icg = 2; irc = 2 ; nu = 1.0;                         % Search direction.
isg_m = 0.05; isg_al0=2; isg_k=0.3;   % stochastic gradient
iheader = 0;
%
%
% Functions for optimization
L   = @(w) norm(y(Xtr,w)-ytr)^2 + (la*norm(w)^2)/2; 
% Optimization
%
t1=clock;
[Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target, tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,ialmax,kmaxBLS,epsal,c1,c2,isd,isg_m,isg_al0,isg_k,icg,irc,nu,iheader);
t2=clock;
fprintf(' wall time = %6.1d s, optimization time = %6.2d s\n', etime(t2,t1),tex);
%


