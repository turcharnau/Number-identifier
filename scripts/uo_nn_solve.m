function [Xtr,ytr,wo,fo,tr_acc,Xte,yte,te_acc,niter,tex]=uo_nn_solve(num_target, tr_freq,tr_seed,tr_p,te_seed,te_q,la,epsG,kmax,ils,almax,kmaxBLS,epsal,c1,c2,isd,isg_m,isg_al0, isg_k, icg,irc,nu, iheader)
%%Generate the training dataset.
[Xtr,ytr] = uo_nn_dataset(tr_seed, tr_p, num_target, tr_freq);
%%Generate the test dataset.
[Xte,yte] = uo_nn_dataset(te_seed, te_q, num_target, tr_freq);


%OPTIMIZATION-------------------------------------------------
t1=clock;
tr_m = floor(isg_m*tr_p);
[wo,fo,niter, dk, alk, iWk,wk]=optimize(Xtr, ytr, la,epsG,kmax,almax,kmaxBLS,epsal,c1,c2,isd,icg,irc,nu, isg_al0, isg_k, tr_m);
t2=clock;
tex = etime(t2,t1);
%------------------------------------------------------------

sig = @(X) 1./(1+exp(-X));
y = @(X,w) sig(w'*sig(X));


sum=0;
for j = 1:tr_p
  sum=sum+(round(y(Xtr(:,j),wo))==ytr(:,j));
end
tr_acc =100/tr_p*sum;




sum=0;
for j = 1:te_q
  sum=sum+(round(y(Xte(:,j),wo))==yte(:,j));
end

te_acc =100/te_q*sum;


% Output

L   = @(w) norm(y(Xtr,w)-ytr)^2 + (la*norm(w)^2)/2; 
gL = @(w) 2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))'+la*w;

fprintf('Training data set generation\n');
fprintf(' numtarget= %2d\n',num_target);
fprintf(' tr_freq= %2d\n',tr_freq);
fprintf(' tr_p= %2d\n',tr_p);
fprintf(' tr_seed= %2d\n',tr_seed);


fprintf('------------------------------------------------\n');

Lk = []; gLk = []; gdk = [];
for k = 1:niter Lk = [Lk,L(wk(:,k))]; gLk=[gLk,gL(wk(:,k))]; end % L(wk) and gL(wk)
for k = 1:niter-1 gdk = [gdk,gLk(:,k)'*dk(:,k)]; end % Descent condition
gdk = [gdk,0];
fprintf('Optimization\n');
fprintf(' L2 reg. lambda = %4',la);
fprintf(' epsG= %3.1e, kmax= %4d\n', epsG, kmax);
fprintf(' almax= %2d, kmaxBLS= %3.1e, epsBLS= %4.2f\n',almax,kmaxBLS,epsal);
fprintf(' c1= %3.2f, c2= %3.2f, isd= %1d\n',c1,c2,isd);
fprintf(" k     al     iW     g'*d     f     ||g||     \n");
for k = 1:2
fprintf('%5d %7.4f %7.4f %3d %+3.1e %4.2e\n', k, alk(k), iWk(k), gdk(k), Lk(k) ,norm(gLk(:,k)));
end
fprintf('------------------------------------------------\n');
for k = 1:niter
fprintf('%5d %7.4f %7.4f %3d %+3.1e %4.2e\n', k, alk(k), iWk(k), gdk(k), Lk(k) ,norm(gLk(:,k)));
end
fprintf(" k     al     iW     g'*d     f     ||g||     \n");
fprintf('------------------------------------------------\n');

disp(['Wo: [' num2str(wo(:).') ']']) ;

fprintf('Test data set generation\n');
fprintf(' te_q= %2d\n',te_q);
fprintf(' te_seed= %2d\n',te_seed);
fprintf(' tr_accuracy= %3.2f\n',tr_acc);
fprintf(' te_accuracy= %3.2f\n',te_acc);


%subplot(2,1,1); 
%uo_nn_Xyplot(Xtr,ytr,wo);
%subplot(2,1,2); 
uo_nn_Xyplot(Xte,yte,wo);

end