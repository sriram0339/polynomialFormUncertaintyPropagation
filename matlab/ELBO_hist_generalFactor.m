function [elbo] = ELBO_hist_generalFactor(theta,f2,w,maxw,minw,n_int,Clusters,indexTab)
    size_C=size(Clusters,2);
    count_c=0;
    for i=1:size_C
        d_C=size(Clusters{i},2);
        num_Cell = n_int^d_C;
        A{i} = theta(count_c+1:count_c+num_Cell).';
        count_c=count_c+num_Cell;
    end
    
%A=reshape(theta,n_int,d).';
% Calculate objective f
%integral over parameter space
syms F
for i=1:size_C
    %delta=(maxw(i)-minw(i))/2/n_int;
    delta=(maxw(Clusters{i})-minw(Clusters{i}))/2/n_int;
    d_C=size(Clusters{i},2);
    num_Cell=n_int^d_C;
    F=0;
    
    for j=1:num_Cell
        
        ind_temp=indexTab{i}(j,:).';

        midw_temp=minw(Clusters{i})+ 2*delta.*(ind_temp-1/2);
        midw_temp=midw_temp.';
        %F = F + int(f2*A(i,j)/dw,w(i),[minw_temp maxw_temp]);
        F = F + subs(f2*A{i}(1,j),w(Clusters{i}),midw_temp);%instead of integrating, we use the center point to represnt the expectation
    end
    f2=F;
end
Hq=0;
for i=1:size_C
    d_C=size(Clusters{i},2);
    %dw=(maxw(i)-minw(i))/n_int;
    Hqi=0;
    for j=1:n_int^d_C
        Hqi=Hqi+A{i}(1,j)*log(A{i}(1,j));%A(i,j)*log(A(i,j));
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