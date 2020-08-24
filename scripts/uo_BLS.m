function [alphas,iout] = uo_BLS(L,g_BLS,wk,dk,almax,c1,c2,kBLSmax,epsal,k,kmax,isg_al0,isg_k,isd,alk)

if (isd==7)
    al_sg = 0.01*isg_al0;

    k_sg = floor(isg_k*kmax);

    if(k<=k_sg) alphas=(1-k/k_sg)*isg_al0+(k/k_sg)*al_sg;
    else alphas = al_sg;
    end
    iout=0;
else
    if k==1, ialmax=1;
    else
        if almax==1, ialmax=alk(k-1)*g_BLS(wk(:,k-1))'*dk(:,k-1)/g_BLS(wk(:,k))'*dk(:,k);
        elseif almax==2, ialmax= 2*(L(wk(:,k))-L(wk(:,k-1)))/(g_BLS(wk(:,k))'*dk(:,k));
        end
    end
    
    [alphas, iout]=uo_BLSNW32(L,g_BLS,wk(:,k),dk(:,k),ialmax,c1,c2,kBLSmax,epsal);
end
end 