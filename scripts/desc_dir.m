function [d] = desc_dir(wk, gk, Bk, isd, icg, irc, nu, dk, k, tr_m)

n = length(wk(:,1));

%gradient method or first iteration of cgm
if isd==1 || (isd==2 && k == 1)

    d = -gk(:,k);

%conjugate gradient method
elseif isd == 2
  %restart condition
  
    if irc == 1 && mod(k,n) == 0 || irc == 2 && ((abs(gk(:,k)'*gk(:,k-1)))/(norm(gk(:,k)))^2)>=nu 

      b = 0;

  %no restart condition nor first iteration
    else

        if icg==1 %Fletcher Reedes
  
            b=(gk(:,k)'*gk(:,k))/(norm(gk(:,k-1)))^2;
        
        elseif icg ==2 %Polak Ribiere 
    
            b=max(0,(gk(:,k)'*(gk(:,k)-gk(:,k-1)))/(norm(gk(:,k-1))^2));
            
        end
    end
       
    d =-gk(:,k)+b.*dk(:,k-1);    

elseif isd == 3
  d=-Bk*gk(:,k);

elseif isd == 7

  d=-(1/tr_m)*gk(:,k);

end
end
