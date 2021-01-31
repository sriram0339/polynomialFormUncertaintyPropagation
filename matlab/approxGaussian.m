%clear all;
addpath('./ModelParamerter')

%parameterSetting_neuronmodel
%parameterSetting_p53
%parameterSetting_ebola
%parameterSetting_honeybee
%parameterSetting_LaubLoomis

tic
logP=readLogP('../outputs/logPoly_Trunc.txt');
logp=0;
for i=1:size(logP,2)
logp = logp + str2sym(logP{i});
end

%logp=45.3331*w3 -5206.07*w3^2 -40.9975*w2 + 2299.16 *w2* w3 -1203.14*w2^2;
% pI = [-9.33038,-9.12161];
% %pI = [-11.0634,-10.884];%p53
% pI = pI - mean(pI); %error;

%find H such that Q = -1/2*[x1,x2]*K*[x1;x2] + h*[x1;x2]
Q=logp;
X=w;
% d=2;
% syms w2 w3;
% Q = 45.3331*w3 -5206.07*w3^2 -40.9975*w2 + 2299.16 *w2* w3 -1203.14*w2^2;
% X = [w2 w3];
% syms w6 w7;
% Q =  -2452.05*w7 -766819*w7^2 -3958.78*w6 -1.6587e+06*w6* w7 -1.16768e+06*w6^2;
% X = [w6 w7];

X0=zeros(1,d);
% h(i) = dQ/dwi|(w=0)
for i=1:d
h(1,i) = subs(diff(Q,X(i)),X,X0); 
end
K = -2*(hessian(Q)/2);

%%%test Q=-1/2*[x1,x2]*K*[x1;x2]+ h*[x1;x2]
%Qnew = expand(-1/2*X*K*X.'+ h*X.')

%convert canonical form into Gaussian form by K=\Sigma^-1, h =\Sigma^-1 * \mu 
K=vpa(K);%using decimals instead of rational fractions
h=vpa(h);
Sigma = inv(K);
mu = Sigma * h';
mu1 = mu + (mintheta+maxtheta)/2;% Q is computed with parameters centered at 0
Qg = -1/2*((X-mu1')*K*(X-mu1').');%Gaussian dist: 1/N*exp^(Qg);

eig_S=eig(double(Sigma));
is_PosSemi=all(eig_S(:)>=0);
if(is_PosSemi)
%plot Gaussian martingale 
mu_=double(mu1');
for i=1:d
    subplot(d,1,i);
    std_(i)=sqrt(double(Sigma(i,i)));
    x = mintheta(i):std_(i)/10:maxtheta(i);
    y = pdf('Normal',x,mu_(i),std_(i));
    plot(x,y)
    hold on
end
end
toc
% % %%compute the coefficients of histogram
% minw=-(maxtheta-mintheta)/2;%-10;
% maxw=(maxtheta-mintheta)/2;%10;
% mu0 = Sigma * h';
% coef=zeros(d,numInt);
% for i=1:d
%     delta=(maxw(i)-minw(i))/2/numInt;
%     subplot(d,1,i);
%     std_=sqrt(double(Sigma(i,i)));
%     x = minw(i)+delta:2*delta:maxw(i)-delta;
%     y = pdf('Normal',x,mu0(i),std_);
%     coef(i,:)=y;
% %     %x=x+(mintheta(i)+maxtheta(i))/2;
% %     plot(x,y)
% %     hold on
% end

% %%compute the coefficients of Taylor series of Gaussian
% D_poly=50;%50;
% coef=zeros(2,2*D_poly+1);
% for i=1:d
%     std_=sqrt(double(Sigma(i,i)));
%     fact=(1/std_/sqrt(2*pi));
%     %exp(x^2)=sum (x^2n)/n!
%     
%     fact_temp=1;
%     coef(i,1)=fact;
%     for j=1:D_poly 
%         fact_temp=fact_temp*j;
%         coef(i,2*j+1)=fact*(-1/2/(std_^2))^j/fact_temp;
%     end    
% end
% mu0 = Sigma * h';
%%plot Taylor series of Gaussian
% for i=1:d
%     subplot(d,1,i);
%     std_=sqrt(double(Sigma(i,i)));
%     x = -2*std_:std_/10:2*std_;
%     %x = mintheta(i):std_/10:maxtheta(i);
%     pow_x=ones(1,size(x,2));
%     y=coef(i,1)*ones(1,size(x,2));
%     for j=1:D_poly 
%         %pow_x(i)
%         pow_x=pow_x.*x.^2;
%         y=y+coef(i,2*j+1).*pow_x;
%     end
%     x=x+mu_(i);
%     plot(x,y)
%     hold on
% end

% %compute the normalization constant
% const_N = (2*3.1415926)^(d/2)*det(Sigma)^(1/2);
% %const_N = (2*3.1415926)^(d/2)*norm(Sigma)^(1/2);
% %compute the prob. on the grid points
% filename = 'centra list.txt';
% [centra,delimiterOut]=importdata(filename);
% filename = 'deltap list.txt';
% [deltap,delimiterOut]=importdata(filename);
% ncell=size(centra,1);
% 
% %offset= (mintheta+maxtheta)/2;%orignal polynomial Q is computed with parameters centered at 0
% logConst_N=log(const_N);
% 
% for i=1:ncell 
%     X_temp = centra(i,:);%centra(i,:) - offset';
%     prob2(i,1) = vpa(subs(Qg,X,X_temp)) - logConst_N ;%+ pI(1);%lower
%     prob2(i,2) = vpa(subs(Qg,X,X_temp)) - logConst_N ;%+ pI(2);%upper
% end
% 
% 
% %project 2d to 1d
% ncell=size(prob2,1);
% den=zeros(ncell,2);
% for i=1:ncell
%     den(i,1)=exp(prob2(i,1))*deltap(i,2)*(2^1);%/deltap(i,1)/(2^1);%lower bound
%     den(i,2)=exp(prob2(i,2))*deltap(i,2)*(2^1);%/deltap(i,1)/(2^1);%upper bound
% end
% 
% centra_1=[centra(1,1)];
% den_1= [den(1,1) den(1,2)];
% deltap_1=[deltap(1,1)];
% for i=2:ncell
%     ind_temp1=find(centra_1(:,1)==centra(i,1));
%     if size(ind_temp1,1)==0
%         
%         centra_1=[centra_1;centra(i,1)];
%         den_1=[den_1; den(i,1) den(i,2)];
%         deltap_1=[deltap_1; deltap(i,1)];
%     else           
%         ind_den_1=ind_temp1;
%         den_1(ind_den_1,:)=den_1(ind_den_1,:)+[den(i,1) den(i,2)];
%        
%     end
% 
% end
% 
% deInd=[];
% for i=1:size(centra_1(:,1),1)
%     for j=1:size(centra_1(:,1),1)
%         if (i~=j)
%             if(centra_1(i,1)+deltap_1(i,1) <=centra_1(j,1)+deltap_1(j,1) && centra_1(i,1)-deltap_1(i,1) >= centra_1(j,1)-deltap_1(j,1))
%                 den_1(i,:)=den_1(i,:)+den_1(j,:);
%                 deInd=[deInd;j];
%             end
%             
%         end
%     end
% end
% 
% 
% %plot(centra_1(:,1),den_1(:,1),'.b',centra_1(:,1),den_1(:,2),'.r')
% deInd_uniq = unique(deInd);
% centra_1new=centra_1(setdiff(1:end,deInd_uniq),1);
% den_1lo=den_1(setdiff(1:end,deInd_uniq),1);
% den_1up=den_1(setdiff(1:end,deInd_uniq),2);
% den_1new=[den_1lo den_1up];
% subplot(2,1,1);
% %plot(centra_1new(:,1),den_1new(:,1),'.b',centra_1new(:,1),den_1new(:,2),'.r')
% nCnew=size(centra_1new,1);
% Sden_1new=zeros(nCnew,2);
% Scentra_1new=zeros(nCnew,1);
% [sortV , sortIndex] = sort(centra_1new(:,1));
% Scentra_1new(:,1)=sortV;
% for i=1:nCnew
%     Sden_1new(i,1) = den_1new( sortIndex(i), 1);
% end
% stairs(Scentra_1new(:,1),Sden_1new(:,1),'*-b','LineWidth',1) 
% hold on;
% for i=1:nCnew
%     Sden_1new(i,2) = den_1new( sortIndex(i), 2);
% end
% stairs(Scentra_1new(:,1),Sden_1new(:,2),'*-r','LineWidth',1)
% 
% 
% ncell=size(prob2,1);
% den=zeros(ncell,2);
% for i=1:ncell
%     den(i,1)=exp(prob2(i,1))*deltap(i,1)*(2^1);%/deltap(i,2)/(2^1);%lower bound
%     den(i,2)=exp(prob2(i,2))*deltap(i,1)*(2^1);%/deltap(i,2)/(2^1);%upper bound
% end
% 
% centra_2=[centra(1,2)];
% den_2= [den(1,1) den(1,2)];
% deltap_2=[deltap(1,2)];
% for i=2:ncell
%     ind_temp1=find(centra_2(:,1)==centra(i,2));
%     if size(ind_temp1,1)==0
%         
%         centra_2=[centra_2;centra(i,2)];
%         den_2=[den_2; den(i,1) den(i,2)];
%         deltap_2=[deltap_2; deltap(i,2)];
%     else           
%         ind_den_2=ind_temp1;
%         den_2(ind_den_2,:)=den_2(ind_den_2,:)+[den(i,1) den(i,2)];
%        
%     end
% 
% end
% 
% deInd=[];
% for i=1:size(centra_2(:,1),1)
%     for j=1:size(centra_2(:,1),1)
%         if (i~=j)
%             if(centra_2(i,1)+deltap_2(i,1) <=centra_2(j,1)+deltap_2(j,1) && centra_2(i,1)-deltap_2(i,1) >= centra_2(j,1)-deltap_2(j,1))
%                 den_2(i,:)=den_2(i,:)+den_2(j,:);
%                 deInd=[deInd;j];
%             end
%             
%         end
%     end
% end
% 
% deInd_uniq = unique(deInd);
% %setdiff(1:end,deInd_uniq)
% centra_2new=centra_2(setdiff(1:end,deInd_uniq),1);
% den_2lo=den_2(setdiff(1:end,deInd_uniq),1);
% den_2up=den_2(setdiff(1:end,deInd_uniq),2);
% den_2new=[den_2lo den_2up];
% 
% 
% subplot(2,1,2);
% %plot(centra_2new(:,1),den_2new(:,1),'.b',centra_2new(:,1),den_2new(:,2),'.r')
% nCnew=size(centra_2new,1);
% Sden_2new=zeros(nCnew,2);
% Scentra_2new=zeros(nCnew,1);
% [sortV , sortIndex] = sort(centra_2new(:,1));
% Scentra_2new(:,1)=sortV;
% for i=1:nCnew
%     Sden_2new(i,1) = den_2new( sortIndex(i), 1);
% end
% stairs(Scentra_2new(:,1),Sden_2new(:,1),'*-b','LineWidth',1) 
% hold on;
% for i=1:nCnew
%     Sden_2new(i,2) = den_2new( sortIndex(i), 2);
% end
% stairs(Scentra_2new(:,1),Sden_2new(:,2),'*-r','LineWidth',1)






% tol_proj=1e-20;
% centra_mean_OP=zeros(numInt,d);
% for k=1:d
%     temp_delta= (maxtheta(k)-mintheta(k))/numInt/2;
%     for i=1:numInt
%         centra_mean_OP(i,k)=mintheta(k) + temp_delta*(2*(i-1)+1);
%     end
% end
% 
% norm_meanDen_low=zeros(numInt,d);
% for k=1:d
%    temp_delta= (maxtheta(k)-mintheta(k))/numInt/2; 
%    for i=1:ncell
%       for j=1:numInt
%          if(abs(centra_mean_OP(j,k)-centra(i,k)) < tol_proj)
%             tempInd = j;
%             break
%          end
%       end
%       norm_meanDen_low(tempInd,k) = norm_meanDen_low(tempInd,k) + exp(prob2(i,1))/2/temp_delta;
%    end
%     
% end
% 
% subplot(2,1,1);
% stairs(centra_mean_OP(:,1),norm_meanDen_low(:,1),'*-b','LineWidth',1) 
% subplot(2,1,2);
% stairs(centra_mean_OP(:,2),norm_meanDen_low(:,2),'*-b','LineWidth',1) 
% 
% hold on
% norm_meanDen_up=zeros(numInt,d);
% for k=1:d
%    temp_delta= (maxtheta(k)-mintheta(k))/numInt/2; 
%    for i=1:ncell
%       for j=1:numInt
%          if(abs(centra_mean_OP(j,k)-centra(i,k)) < tol_proj)
%             tempInd = j;
%             break
%          end
%       end
%       norm_meanDen_up(tempInd,k) = norm_meanDen_up(tempInd,k) + exp(prob2(i,2))/2/temp_delta;
%    end
%     
% end
% 
% subplot(2,1,1);
% stairs(centra_mean_OP(:,1),norm_meanDen_up(:,1),'*-r','LineWidth',1) 
% subplot(2,1,2);
% stairs(centra_mean_OP(:,2),norm_meanDen_up(:,2),'*-r','LineWidth',1) 
% 
% 
% for i=1:d
%    temp_delta= (maxtheta(i)-mintheta(i))/numInt/2;
% 
%  Ep=sum(centra_mean_OP(:,i)'*norm_meanDen_up(:,i)*temp_delta*2);
%  Ep
% end
%  




% filename = 'theta_sample.txt';
% [theta1,delimiterOut]=importdata(filename);
% N=size(theta1,1);
% % DISPLAY SAMPLING DYNAMICS
% theta1=theta1';
% x1=theta1';
% subplot(2,1,1);
% [n2 x2] = hist(x1(:,1), ceil(sqrt(N))); 
% bar(x2, n2/(N*(x2(2)-x2(1))));   hold on;               % Normalized histogram
% xlabel('p(alpha)', 'FontSize', 15) 
% subplot(2,1,2);
% [n2 x2] = hist(x1(:,2), ceil(sqrt(N))); 
% bar(x2, n2/(N*(x2(2)-x2(1))));   hold on;               % Normalized histogram 
% xlabel('p(beta)', 'FontSize', 15)

% mu=double(mu');
% Sigma=double(Sigma);
% x = mintheta(1):(maxtheta(1)-mintheta(1))/50:maxtheta(1); %// x axis
% y = mintheta(2):(maxtheta(2)-mintheta(2))/50:maxtheta(2); %// y axis
% % x = -5*theta(d+1):theta(d+1)/10:5*theta(d+1); %// x axis
% % x = mu(1)+x;
% % y = -5*theta(d+2):theta(d+2)/10:5*theta(d+2); %// y axis
% % y = mu(2)+y;
% [X Y] = meshgrid(x,y);
% Z = mvnpdf([X(:) Y(:)],mu,Sigma); %// compute Gaussian pdf
% Z = reshape(Z,size(X));
% contour(X,Y,Z), axis equal  %// contour plot; set same scale for x and y...
% %contour(X,Y,Z,'ShowText','on')
% %surf(X,Y,Z,'LineStyle','none') %// ...
