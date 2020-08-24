function [wo,fo,niter, dk, alk, iout,wk]=optimize(Xtr, ytr,la,epsG,kmax,almax,kBLSmax,epsal,c1,c2,isd,icg,irc,nu, isg_al0, isg_k, tr_m)

sig = @(X) 1./(1+exp(-X));
y = @(X,w) sig(w'*sig(X));
L   = @(w) norm(y(Xtr,w)-ytr)^2 + (la*norm(w)^2)/2; 
gL = @(w,X,Y) 2*sig(X)*((y(X,w)-Y).*y(X,w).*(1-y(X,w)))'+la*w;
gL_BLS = @(w) 2*sig(Xtr)*((y(Xtr,w)-ytr).*y(Xtr,w).*(1-y(Xtr,w)))'+la*w;

w=zeros(35,1);

n=35;alk=zeros(1,kmax); iout=zeros(1, kmax);
Bk=zeros(n,n,kmax);wk=zeros(n,kmax); dk=zeros(n,kmax); gk=zeros(n,kmax);
k=1; wk(:,1)=w; Bk(:, :, 1) = eye(n);

if isd~=7
    gk(:,k)=gL(wk(:,k),Xtr,ytr);
    tr_m=1;
else 
    p=size(Xtr,2);
    batch = randi([1,p],1,tr_m);    
    Xs=Xtr(:,batch); ys=ytr(:,batch); 
    gk(:,k)=gL(wk(:,k),Xs,ys);
end

while(k<kmax && norm(gk(:,k))>epsG)
    
  [dk(:,k)] = desc_dir(wk, gk, Bk(:, :, k), isd, icg, irc, nu, dk, k,tr_m);
    
  [alk(:,k),iout(:,k)]=uo_BLS(L,gL_BLS,wk,dk,almax,c1,c2,kBLSmax,epsal,k,kmax,isg_al0,isg_k,isd,alk);
  
  wk(:,k+1) = wk(:,k)+alk(:,k)*dk(:,k);

  if isd~=7
      gk(:,k+1)=gL(wk(:,k+1),Xtr,ytr);
  else 
      p=size(Xtr,2);
      batch = randi([1,p],1,tr_m);    
      Xs=Xtr(:,batch); ys=ytr(:,batch); 
      gk(:,k+1)=gL(wk(:,k+1),Xs,ys);
  end

  if(isd==3)  
    s=wk(:,k+1)-wk(:,k); yk=gk(:,k+1)-gk(:,k); p=1/(yk'*s);
    Bk(:, :, k+1)=(eye(n)-p*s*yk')*Bk(:,:,k)*(eye(n)-p*yk*s')+p*s*s';
  end
  
    k=k+1;
    
end

wo = wk(:,k);
fo = L(wo);
niter = k;

end
