clear all
addpath('./ModelParamerter')

parameterSetting_neuronmodel
%parameterSetting_p53
%parameterSetting_ebola
%parameterSetting_honeybee
%parameterSetting_LaubLoomis
logP=readLogP('../outputs/logPoly.txt');

v = sym('v',[d 1]);

logp=0;
for i=1:size(logP,2)
logp = logp + str2sym(logP{i});
end
%logp = 45.3331*w3 -5206.07*w3^2 -40.9975*w2 + 2299.16 *w2* w3 -1203.14*w2^2;

% mu=sym('mu', [1 d]);
% sigma = sym('sigma', [1 d]);
A = sym('A', [d d]);
B = sym('B', [d 1]);
sqrt2pi_inv = 0.39894;%1/ sqrt(2*pi);
f =logp;
%f =  -2452.05*w7 -766819*w7^2 -3958.78*w6 -1.6587e+06*w6* w7 -1.16768e+06*w6^2;
% string_ = fileread('logPoly.txt');
% count=0;
% for i=1:strlength(string_)
%     if(strcmp(string_(i),'+')||strcmp(string_(i),'-'))
%         count=count+1;
%     end
% end
f = subs(f,w.',(A*v+B));%replace parameters with new variables by a linear transformation
for i=1:d
   
   f= f *  sqrt2pi_inv * exp(-v(i)^2/ 2); 
end

if(is_symbolic==1)
fprintf('compute the integral of f \n')
tic
%integral over parameter space
for i=1:d
    i
    %f = int(f,v(i),[minw(i) maxw(i)]);
    f = int(f,v(i),1000*[minw(i) maxw(i)]);
    %f = int(f,v(i),[-10 10]);
end
toc

%Convert symbolic expression to function handle or file
Anew=reshape(A.',d*d,1);
Anew=[Anew;B];
%matlabFunction(f,'File','expectation_q','Vars',{Anew});
expectation_q=matlabFunction(f,'Vars',{Anew});

%Gradients of expectation of log likelihood f
if(is_grad==1)
 gradF2g = gradient(f,Anew);
 gradF2=matlabFunction(gradF2g,'Vars',{Anew});
%matlabFunction(gradF2g,'File','gradF2','Vars',{Anew});
end

%entropy of Q
Hq = log(det(A));
%matlabFunction(-Hq,'File','negativeH','Vars',{Anew});
negativeH = matlabFunction(-Hq,'Vars',{Anew});
if(is_grad==1)
gradHq = gradient(-Hq,Anew);
gradNegativeH = matlabFunction(gradHq,'Vars',{Anew});
%gradNegativeH = matlabFunction(gradHq,'Vars',{Anew});
end

else
    Anew=reshape(A.',d*d,1);
    Anew=[Anew;B];
end%is_symbolic==1

if(is_symbolic==1)
    if(is_grad==0)
    %sym_theta=[A,B];
    fun = @(theta)  double(-vpa(subs(f,Anew,theta))-vpa(subs(Hq,Anew,theta)));
    %fun = @(theta) -expectation_q(theta) + negativeH(theta);
    else
    fun =@(theta) ELBOwithgrad(theta,expectation_q,negativeH,gradF2,gradNegativeH);
    end%if(is_grad==0)
else
    fun =@(theta) ELBO_multG(theta,f,Anew,v,maxw,minw,n_int_MG,d);
end


%A0 = [0.001 0 ; 0 0.00167];
%A0 = [sqrt(0.000526) 0 ; 0 sqrt(0.0001217)];
% A0 = [0.001 0 ; 0 0.001];
% A_= zeros(4*d,d*d+d);
% b_=zeros(2*d,1);
% for i=1:d
%     A_(i,(i-1)*d+i)=1;
%     A_(d+i,(i-1)*d+i)=-1; 
%     A_(2*d+i,d*d+i)=1;
%     A_(3*d+i,d*d+i)=-1;
%     b_(i)=(maxw(i)-minw(i))/5;
%     b_(d+i)=0;
%     b_(2*d+i)=(maxw(i)-minw(i))/2;
%     b_(3*d+i)=-(maxw(i)-minw(i))/2;
% end

%%% use the information of range of parameter to bound the A,b matrix
lb=zeros(d*d+d,1);
ub=zeros(d*d+d,1);

Alb=zeros(d,d);
Aub=zeros(d,d);
for i=1:d
    Alb(i,:)=-(maxw(i)-minw(i))/5;
    Alb(i,i)=10^-10;
    Aub(i,:)=(maxw(i)-minw(i))/5;
end

lb(1:d*d,1)=reshape(Alb.',d*d,1);
ub(1:d*d,1)=reshape(Aub.',d*d,1);
lb(d*d+1:end,1)=-(maxw-minw)/2;
ub(d*d+1:end,1)=(maxw-minw)/2;

approxGaussian

minfval=10^10;
for k=1:1
  A0 = zeros(d,d);  
  B0 = zeros(d,1);%[0.00; 0.00];
  if(k==1)
    if(is_PosSemi==0)
        mu_=zeros(1,d);
        std_=(maxw.'-minw.')/20;
        mu_=double(mu');
        for i=1:d
        if(double(Sigma(i,i))>0)
            std_(i)=sqrt(double(Sigma(i,i)));
        end
        end%i
    
    end
    B0 = mu_.';
    for i=1:d
        A0(i,i)= std_(i);
    end
    
  else
    

    for j=1:d
        A0(j,j)= rand*(maxw(j)-minw(j))/10;
    end
  end%k==1

%B0 = [-0.016321; 0.0007499];
% B0 = [-0.0024; 0.001];
A0new=reshape(A0.',d*d,1);
theta0= [A0new;B0];%[A0,B0];
%theta0= [A0,B0];
tic
%[theta,fval] = fminsearch(fun, theta0);
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
%[theta,fval] = fminunc(fun, theta0,options);
%[theta,fval] = fminunc(fun, theta0);
if(is_grad==0)
[theta,fval] = fmincon(fun,theta0,[],[],[],[],lb,ub);    
else
[theta,fval] = fmincon(fun,theta0,[],[],[],[],lb,ub,[],options);
end
toc
fval
    if(fval<minfval)
        minfval=fval;
        opt_theta=theta;
    end
end
theta=opt_theta;
% A_opt = theta(:,1:d);
% B_opt = theta(:,d+1:end);
A_opt = reshape(theta(1:d*d,1),d,d).';
B_opt = theta(d*d+1:end,1);
%plot Gaussian martingale 
mu_= B_opt;
Sigma = A_opt*A_opt';
for i=1:d
    subplot(d,1,i);
    std_=sqrt(double(Sigma(i,i)));
    x = -5*std_:std_/10:5*std_;
    x = mu_(i)+x;
    y = pdf('Normal',x,mu_(i),std_);
    x=x+(mintheta(i)+maxtheta(i))/2;
    plot(x,y)
    hold on
end
save opt_theta_multiG opt_theta
% %%% plot contour
% mu = B_opt'; %// data
% Sigma =  A_opt*A_opt'; %// data
% std_=sqrt(double(Sigma(1,1)));
% x = -5*std_:std_/10:5*std_; %// x axis
% x = mu(1)+x;
% std_=sqrt(double(Sigma(2,2)));
% y = -5*std_:std_/10:5*std_; %// y axis
% y = mu(2)+y;
% [X Y] = meshgrid(x,y);
% % F = @(qx,qy) pdf('Normal',qx,theta(1),theta(d+1)).*pdf('Normal',qy,theta(2),theta(d+2));
% % qz = F(X,Y);
% % surf(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,qz,'LineStyle','none');
% Z = mvnpdf([X(:) Y(:)],mu,Sigma); %// compute Gaussian pdf
% Z = reshape(Z,size(X));
% contour(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z), axis equal  %// contour plot; set same scale for x and y...
% %contour(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z,'ShowText','on')
% %surf(X+(mintheta(1)+maxtheta(1))/2,Y+(mintheta(2)+maxtheta(2))/2,Z,'LineStyle','none') %// ... or 3D plot
% 
