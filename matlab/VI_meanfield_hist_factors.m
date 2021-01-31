clear all
%approxGaussian

addpath('./ModelParamerter')

%parameterSetting_neuronmodel
%parameterSetting_p53
parameterSetting_ebola
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
for i=1:size_C
    d_C=size(Clusters{i},2);
    A{i}=sym(char(96+i), [1 n_int^d_C]);
end
sqrt2pi_inv = 0.39894;%1/ sqrt(2*pi);
f =logp;


%% indexTab: generate all combinations for the cell in each Cluster 
for i=1:size_C
 clear L   
 for j=1:size(Clusters{i},2)
     L{j} = [1:n_int];
 end
n = length(L);
[L{:}] = ndgrid(L{end:-1:1});
L = cat(n+1,L{:});
L = fliplr(reshape(L,[],n));
indexTab{i} = L;
end

Anew=[];
for i=1:size_C
    Anew=[Anew;A{i}.'];%reshape(A.',d*n_int,1);
end

if(is_symbolic==1)
fprintf('compute the integral of f \n')
tic
%integral over parameter space
f2=f;
syms F
for i=1:size_C
    i
    delta=(maxw(Clusters{i})-minw(Clusters{i}))/2/n_int;
%     dw=2*delta;
    d_C=size(Clusters{i},2);
    F=0;
    for j=1:n_int^d_C
        ind_temp=indexTab{i}(j,:).';
        
        midw_temp=minw(Clusters{i})+ 2*delta.*(ind_temp-1/2);
        midw_temp=midw_temp.';
%         minw_temp=minw(Clusters{i})+ 2*delta*(ind_temp-1);
%         maxw_temp=minw(Clusters{i})+ 2*delta*(ind_temp);
        %F = F + int(f2*A(i,j)/dw,w(i),[minw_temp maxw_temp]);
        F = F + subs(f2*A{i}(1,j),w(Clusters{i}),midw_temp);
        %F = F + subs(f2*A(i,j),w(i),(minw_temp+maxw_temp)/2);%instead of integrating, we use the center point to represnt the expectation
    end
    f2=F;
end
toc
%Convert symbolic expression to function handle or file

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
for i=1:size_C
    d_C=size(Clusters{i},2);
    %dw=(maxw(i)-minw(i))/n_int;
    Hqi=0;
    for j=1:n_int^d_C
        Hqi=Hqi+A{i}(1,j)*log(A{i}(1,j));%A(i,j)*log(A(i,j));
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
count_c=0;
for i=1:size_C
    d_C=size(Clusters{i},2);
    num_Cell = n_int^d_C;
    for j=1:num_Cell
        gradHq(count_c+j,1)=1+log(A{i}(1,j));%gradHq((i-1)*n_int+j,1)=1+log(A(i,j));
    end
    count_c=count_c+num_Cell;
end
%matlabFunction(gradHq,'File','gradNegativeH','Vars',{Anew});
gradNegativeH = matlabFunction(gradHq,'Vars',{Anew});
end 

if(is_symbolic==1)
    %sym_theta= A;
    if(is_grad==1)
    fun =@(theta) ELBOwithgrad(theta,expectation_q,negativeH,gradF2,gradNegativeH);
    else
    %fun =@(theta) ELBOwithgrad(theta);
    fun = @(theta) -expectation_q(theta) + negativeH(theta);
    end
else
    fun = @(theta) ELBO_hist_generalFactor(theta,f,w,maxw,minw,n_int,Clusters,indexTab);
end

% Aeq_temp=zeros(d,n_int);
% for i=1:size_C
%     d_C=size(Clusters{i},2);
%     num_Cell = n_int^d_C;
%     for j=1:num_Cell
%         Aeq_temp(i,j)=1;
%     end
% end

%total variational variables
num_vc=size(Anew,1);
Aeq=zeros(size_C,num_vc);
count_c=0;
for i=1:size_C
    d_C=size(Clusters{i},2);
    num_Cell = n_int^d_C;
    Aeq_temp=ones(1,num_Cell);
    Aeq(i,count_c+1:count_c+num_Cell)=Aeq_temp(1,:);
    %Aeq(i,(i-1)*n_int+1:i*n_int)
    count_c=count_c+num_Cell;
end
beq=ones(size_C,1);
A_=[Aeq;-Aeq];
db=0.05*ones(size_C,1);
b_=[beq+db;-(beq-db)];
lb = zeros(num_vc,1);
ub = ones(num_vc,1);
%theta0= 1/n_int*ones(d,n_int);%coef;

minfval=10^10;
for k=1:1%10%20
    k
%if(k==1)
 %theta0 = rand(d,n_int);
%reshape theta0 into a vector [a00 a01 ...a0k a10 a11 ...ak1]
theta0=rand(num_vc,1);%reshape(theta0',d*n_int,1);
% else
%     theta0=theta;
% end
% load('opt_theta.mat');
% theta0=opt_theta;
tic


scaleFact=Aeq*theta0;
    count_c=0;
    for j=1:size_C
        num_Cell = size(A{j},2);
        theta0(count_c+1:count_c+num_Cell) = theta0(count_c+1:count_c+num_Cell)./scaleFact(j);
        count_c=count_c+num_Cell;
    end
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
    count_c=0;
    for j=1:size_C
        num_Cell = size(A{j},2);
        theta(count_c+1:count_c+num_Cell) = theta(count_c+1:count_c+num_Cell)./scaleFact(j);
        count_c=count_c+num_Cell;
    end
% scaleFact=Aeq*theta;
% theta1=reshape(theta,n_int,d).';
% theta1=theta1./scaleFact;


count_c=0;
for i=1:size_C
    y_set=zeros(size(Clusters{i},2),n_int);
    num_Cell = size(A{i},2);
    theta_temp = theta(count_c+1:count_c+num_Cell);
    for j=1:num_Cell
        ind_temp=indexTab{i}(j,:).';
        for k=1:size(Clusters{i},2)
            y_set(k,ind_temp(k)) =  y_set(k,ind_temp(k)) +theta_temp(j);
        end
    end
    
    for j=1:size(Clusters{i},2)
    %Clusters{i}
    delta=(maxw(Clusters{i}(j))-minw(Clusters{i}(j)))/2/n_int;
    subplot(d,1,Clusters{i}(j));
    x = minw(Clusters{i}(j))+delta:2*delta:maxw(Clusters{i}(j))-delta;
    
    y = y_set(j,:)/2/delta;

    x=x+(mintheta(Clusters{i}(j))+maxtheta(Clusters{i}(j)))/2;
%plot(x,y)
%bar(x,y)
%stairs(x-delta,y,'-g','LineWidth',1)%move the stair to the center point
    stairs(x,y,'-g','LineWidth',1)
    hold on
    end
    count_c=count_c+num_Cell;
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
