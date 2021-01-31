function [elbo] = ELBO_hist(theta,f2,Anew,v,maxw,minw,n,d)
f2 = subs(f2,Anew,theta);
A = reshape(theta(1:d*d,1),d,d).';
% Calculate objective f
%integral over parameter space
syms F
for i=1:d
    int_bound=10*[minw(i) maxw(i)];
    delta=(int_bound(2)-int_bound(1))/2/n;
    F=0;
    for k=1:n
    minw_temp=int_bound(1)+ 2*delta*(k-1);
    maxw_temp=int_bound(1)+ 2*delta*(k);
    F = F + subs(f2,v(i),(minw_temp+maxw_temp)/2)*2*delta;%int(f,v(i),1000*[minw(i) maxw(i)]);

    end
    f2=F;
end
Hq = log(det(A));
f2
elbo = -double(vpa(f2)) - Hq;

end

