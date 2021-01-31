clear all
addpath('./ModelParamerter')

parameterSetting_neuronmodel
%parameterSetting_p53
%parameterSetting_ebola
%parameterSetting_honeybee
%parameterSetting_LaubLoomis
logP=readLogP('../outputs/logPoly.txt');
logp=0;
for i=1:size(logP,2)
logp = logp + str2sym(logP{i});
end


%replace parameter variables as a symbolic vector w=[w1 w2 ... wd];
%w = sym('w',[1 d]);
mu=sym('mu', [1 d]);
sigma = sym('sigma', [1 d]);
sqrt2pi_inv = 0.39894;%1/ sqrt(2*pi);
f =logp;
%f =  -2452.05*w7 -766819*w7^2 -3958.78*w6 -1.6587e+06*w6* w7 -1.16768e+06*w6^2;

for i=1:d
   f= f *  sqrt2pi_inv * (1/sigma(i)) * exp(-(w(i)-mu(i))^2/ 2/(sigma(i)^2)); 
   %f = f * (1/sigma(i)/sqrt(2*pi)) * exp(-(w(i)-mu(i))^2/ 2/(sigma(i)^2));
end

fprintf('compute the integral of f \n')
tic
%integral over parameter space
for i=1:d
    f = int(f,w(i),1000*[minw(i) maxw(i)]);
    %f = int(f,w(i),[-10 10]);%int(f,w(i),[-1 1]);
end
toc
%Anew=[mu,sigma].';
%expectation_q=matlabFunction(f,'Vars',{Anew});

% gradF2g = gradient(f,Anew);
% gradF2=matlabFunction(gradF2g,'Vars',{Anew});
%negative entropy of Q
% Hq=0;
% for i=1:d
%    
%    qi = (1/sigma(i)/sqrt(2*pi)) * exp(-(w(i)-mu(i))^2/ 2/(sigma(i)^2))* (log(1/sigma(i)/sqrt(2*pi)) - (w(i)-mu(i))^2/ 2/(sigma(i)^2));
%    Hqi= int(qi,w(i),[-1 1]);
%    Hq=Hq+Hqi;
% end

% qi=1;
% logqi=0;
% for i=1:d
%    qi=qi * (1/sigma(i)/sqrt(2*pi)) * exp(-(w(i)-mu(i))^2/ 2/(sigma(i)^2));
%    logqi = logqi+  (log(1/sigma(i)/sqrt(2*pi)) - (w(i)-mu(i))^2/ 2/(sigma(i)^2));
% end
% qi = logqi*qi;
% Hq=qi;
% for i=1:d
%    Hq= int(Hq,w(i),[-1 1]);
% end
% negativeH =-sum(0.5*log(2*pi*Anew(d+1:end).^2)+0.5);
% negativeH = matlabFunction(negativeH,'Vars',{Anew});
% 
% gradHq = gradient(negativeH,Anew);
% gradNegativeH = matlabFunction(gradHq,'Vars',{Anew});

sym_theta= [mu,sigma].';
%fun = @(theta)  -vpa(subs(f,sym_theta,theta))+vpa(subs(Hq,sym_theta,theta));
%fun = @(theta)  ELBOwithgrad_uG(theta,sym_theta,f,negativeH,gradF2g,gradNegativeH);
fun = @(theta)  -double(vpa(subs(f,sym_theta,theta)))-sum(0.5*log(2*pi*theta(d+1:end).^2)+0.5);
%fun = @(theta)  -expectation_q-sum(0.5*log(2*pi*theta(d+1:end).^2)+0.5);


lb=zeros(1,2*d).';
ub=zeros(1,2*d).';
lb(1:d)=-(maxw-minw)/2;
lb(d+1:2*d)=0;
ub(1:d)=(maxw-minw)/2;
ub(d+1:2*d)=(maxw-minw);%(maxw-minw)/5;
% theta0= [-0.0024,0.001,0.001,0.00167];
approxGaussian
if(is_PosSemi)
    theta0= [double(mu'), std_].';
else
    mu_=zeros(1,d);
    std_=(maxw.'-minw.')/20;
    mu_=double(mu');
    for i=1:d
        if(double(Sigma(i,i))>0)
            std_(i)=sqrt(double(Sigma(i,i)));
        end
    end
    theta0= [mu_,std_].';
end
tic
%options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
options = optimoptions('fmincon','SpecifyObjectiveGradient',false);

nonlcon = [];
%[theta,fval] = fminsearch(fun, theta0);
%[theta,fval] = fminunc(fun, theta0);
%[theta,fval] = fmincon(fun,theta0);
[theta,fval] = fmincon(fun,theta0,[],[],[],[],lb,ub,nonlcon,options);
toc
for i=1:d
subplot(d,1,i);
% x = minw(i):theta(d+i)/10:maxw(i);
x = -theta(d+i)*5:theta(d+i)/10:theta(d+i)*5;
x = theta(i)+x;
y = pdf('Normal',x,theta(i),theta(d+i));
%x=x+0.5;
x=x+(mintheta(i)+maxtheta(i))/2;
plot(x,y)
end


% mu = [theta(1),theta(2)]; %// data
% sigma = [theta(d+1)^2 0; 0 theta(d+2)^2]; %// data
% x = -5*theta(d+1):theta(d+1)/10:5*theta(d+1); %// x axis
% x = theta(1)+x;
% y = -5*theta(d+2):theta(d+2)/10:5*theta(d+2); %// y axis
% y = theta(2)+y;
% [X Y] = meshgrid(x,y);
% % F = @(qx,qy) pdf('Normal',qx,theta(1),theta(d+1)).*pdf('Normal',qy,theta(2),theta(d+2));
% % qz = F(X,Y);
% % surf(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,qz,'LineStyle','none');
% Z = mvnpdf([X(:) Y(:)],mu,sigma); %// compute Gaussian pdf
% Z = reshape(Z,size(X));
% contour(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z), axis equal  %// contour plot; set same scale for x and y...
% %contour(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z,'ShowText','on')
% %surf(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z,'LineStyle','none') %// ... or 3D plot
