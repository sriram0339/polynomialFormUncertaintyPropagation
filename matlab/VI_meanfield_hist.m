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
%logp =str2sym(logP);
% fid = fopen('logPtest.txt','wt');
% fprintf(fid, '%s', logP);
% fclose(fid);

% fileID = fopen('logPtest.txt','r');
% S = fscanf(fileID,'%s')
%replace parameter variables as a symbolic matrix;
A=sym('A', [d n_int]);
sqrt2pi_inv = 0.39894;%1/ sqrt(2*pi);
f =logp;

if(is_symbolic==1)
fprintf('compute the integral of f \n')
tic
%integral over parameter space
f2=f;
syms F
for i=1:d
    i
    delta=(maxw(i)-minw(i))/2/n_int;
    dw=2*delta;
    
    F=0;
    for j=1:n_int
        
        minw_temp=minw(i)+ 2*delta*(j-1);
        maxw_temp=minw(i)+ 2*delta*(j);
        %F = F + int(f2*A(i,j)/dw,w(i),[minw_temp maxw_temp]);
        F = F + subs(f2*A(i,j),w(i),(minw_temp+maxw_temp)/2);%instead of integrating, we use the center point to represnt the expectation
    end
    f2=F;
end
toc
%Convert symbolic expression to function handle or file
Anew=reshape(A.',d*n_int,1);
%matlabFunction(f2,'File','expectation_q','Vars',{Anew});
expectation_q=matlabFunction(f2,'Vars',{Anew});
end%is_symbolic==1
if(is_grad==1)
%Gradients of expectation of log likelihood f2
gradF2g = gradient(f2,Anew);
%matlabFunction(gradF2g,'File','gradF2','Vars',{Anew});
gradF2=matlabFunction(gradF2g,'Vars',{Anew});    
end%is_grad==1

if(is_symbolic==1)
%negative entropy of Q
syms Hq Hqi
Hq=0;
for i=1:d
    dw=(maxw(i)-minw(i))/n_int;
    Hqi=0;
    for j=1:n_int
        Hqi=Hqi+A(i,j)*log(A(i,j));
    end
   Hq= Hq+Hqi;
end
%matlabFunction(Hq,'File','negativeH','Vars',{Anew});
negativeH = matlabFunction(Hq,'Vars',{Anew});
end

if(is_grad==1)
%Gradients of the negative entropy of Q
% syms gradHq
% gradHq=zeors(d*n_int,1);
for i=1:d
    for j=1:n_int
        gradHq((i-1)*n_int+j,1)=1+log(A(i,j));
    end
end
%matlabFunction(gradHq,'File','gradNegativeH','Vars',{Anew});
gradNegativeH = matlabFunction(gradHq,'Vars',{Anew});
end 

if(is_symbolic==1)
    sym_theta= A;
    if(is_grad==1)
    fun =@(theta) ELBOwithgrad(theta,expectation_q,negativeH,gradF2,gradNegativeH);
    else
    %fun =@(theta) ELBOwithgrad(theta);
    fun = @(theta) -expectation_q(theta) + negativeH(theta);
    end
else
    fun = @(theta) ELBO_hist(theta,f,w,maxw,minw,n_int,d);
end

Aeq_temp=zeros(d,n_int);
for i=1:d
    for j=1:n_int
        Aeq_temp(i,j)=1;
    end
end
Aeq=zeros(d,d*n_int);
for i=1:d
    Aeq(i,(i-1)*n_int+1:i*n_int)=Aeq_temp(i,:);
end
beq=ones(d,1);
A_=[Aeq;-Aeq];
db=0.05*ones(d,1);
b_=[beq+db;-(beq-db)];
lb = zeros(d*n_int,1);
ub = ones(d*n_int,1);
%theta0= 1/n_int*ones(d,n_int);%coef;

minfval=10^10;
for k=1:2%10%20
    k
%if(k==1)
 theta0 = rand(d,n_int);
%reshape theta0 into a vector [a00 a01 ...a0k a10 a11 ...ak1]
theta0=reshape(theta0',d*n_int,1);
%  else
%      theta0=theta;
%  end
% load('opt_theta.mat');
% theta0=opt_theta;
tic


scaleFact=Aeq*theta0;
theta0=reshape(theta0,n_int,d).';
theta0=theta0./scaleFact;
theta0=reshape(theta0',d*n_int,1);
Aeq*theta0;
opts=optimoptions('fmincon', 'MaxFunctionEvaluations',10000);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',10000);
nonlcon = [];
%[theta,fval] = fmincon(fun,theta0,[],[],Aeq,beq,lb,ub);

%[theta,fval] = fmincon(fun,theta0,[],[],Aeq,beq,lb,ub,nonlcon,options);
if(is_grad==1)
[theta,fval] = fmincon(fun,theta0,A_,b_,[],[],lb,ub,nonlcon,options);%out of memory for honeybee
else
[theta,fval] = fmincon(fun,theta0,A_,b_,[],[],lb,ub,[],opts);
end
toc
fval
    if(fval<minfval)
        minfval=fval;
        opt_theta=theta;
    end
end%k
theta=opt_theta;
scaleFact=Aeq*theta;
theta1=reshape(theta,n_int,d).';
theta1=theta1./scaleFact;

for i=1:d
    delta=(maxw(i)-minw(i))/2/n_int;
subplot(d,1,i);
x = minw(i)+delta:2*delta:maxw(i)-delta;
y = theta1(i,:)/2/delta;

x=x+(mintheta(i)+maxtheta(i))/2;
%plot(x,y)
%bar(x,y)
%stairs(x-delta,y,'-g','LineWidth',1)%move the stair to the center point
stairs(x,y,'-g','LineWidth',1)
hold on

end
save opt_theta opt_theta;
% for i=1:6
%     subplot(3,2,i);
%   delta=(maxw(i)-minw(i))/2/n_int;
% x = minw(i)+delta:2*delta:maxw(i)-delta;
% y = theta1(i,:)/2/delta;
% 
% x=x+(mintheta(i)+maxtheta(i))/2;
% stairs(x,y,'-g','LineWidth',1)
% hold on
% 
% end
% 
% for i=7:d
% subplot(4,2,i-6);
%   delta=(maxw(i)-minw(i))/2/n_int;
% x = minw(i)+delta:2*delta:maxw(i)-delta;
% y = theta1(i,:)/2/delta;
% 
% x=x+(mintheta(i)+maxtheta(i))/2;
% 
% stairs(x,y,'-g','LineWidth',1)
% hold on
% 
% end

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
