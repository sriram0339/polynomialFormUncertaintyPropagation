function [elbo] = ELBO_hist(theta,f2,w,maxw,minw,n_int,d)
A=reshape(theta,n_int,d).';
% Calculate objective f
%integral over parameter space
syms F
for i=1:d
    delta=(maxw(i)-minw(i))/2/n_int;
%     dw=2*delta;
    
    F=0;
    for j=1:n_int
        
        minw_temp=minw(i)+ 2*delta*(j-1);
        maxw_temp=minw(i)+ 2*delta*(j);
        %F = F + int(f2*A(i,j)/dw,w(i),[minw_temp maxw_temp]);
        F = F + subs(f2*A(i,j),w(i),(minw_temp+maxw_temp)/2);%instead of integrating, we use the center point to represnt the expectation
    end
    f2=F;
end
Hq=0;
for i=1:d
    Hqi=0;
    for j=1:n_int
        Hqi=Hqi+A(i,j)*log(A(i,j));
    end
   Hq= Hq+Hqi;
end
elbo = -double(vpa(f2)) + Hq;

end
%(f2(w1_1,w2_1)*A(1,1)+f2(w1_2,w2_1)*A(1,2))*A(2,1)+(f2(w1_1,w2_2)*A(1,1)+f2(w1_2,w2_2)*A(1,2))*A(2,2)
%/d A(1,1)
%f2(w1_1,w2_1)*A(2,1)+f2(w1_1,w2_2)*A(2,2)
%/d A(1,2)
%f2(w1_2,w2_1)*A(2,1)+f2(w1_2,w2_2)*A(2,2)
%/d A(2,1)
%(f2(w1_1,w2_1)*A(1,1)+f2(w1_2,w2_1)*A(1,2))
%/d A(2,2)
%(f2(w1_1,w2_2)*A(1,1)+f2(w1_2,w2_2)*A(1,2))
%/d A(1,j)
%f2(w1_j,w2_1)*A(2,1)+f2(w1_j,w2_2)*A(2,2)
%/d A(2,j)
%f2(w1_1,w2_j)*A(1,1)+f2(w1_2,w2_j)*A(1,2)
%/d A(i,j)
%wi:=w_indj, sum_k f2(wi=w_indj,w) 