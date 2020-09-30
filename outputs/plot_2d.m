d=2;%dimension
filename = 'centra list.txt';
[centra,delimiterOut]=importdata(filename);
filename = 'deltap list.txt';
[deltap,delimiterOut]=importdata(filename);

% filename = 'prob list.txt';
% [prob,delimiterOut]=importdata(filename);
filename = 'prob list2.txt';
[prob2,delimiterOut]=importdata(filename);

ncell=size(prob2,1);
den=zeros(ncell,2);
for i=1:ncell
    den(i,1)=exp(prob2(i,1))/deltap(i,1)/deltap(i,2)/(2^d);%lower bound
    den(i,2)=exp(prob2(i,2))/deltap(i,1)/deltap(i,2)/(2^d);%upper bound
end

%%Exp
for j=1:d
    sum1=0;
for i=1:ncell
    
    sum1=sum1+exp(prob2(i,2))*centra(i,j);%upper bound
end
Exp=sum1/sum(exp(prob2(:,2)))
end
%find maximum cell
% [value, ind]=max(den(:,1));
% centra(ind,:)
% [value, ind]=max(den(:,2));
% centra(ind,:)


%plot3(centra(:,1),centra(:,2),den(:,1),'.',centra(:,1),centra(:,2),den(:,2),'.')

%set(gca, 'ZScale', 'log')

for i=1:ncell
    if ((prob2(i,1))>(prob2(i,2)))
        i

    end
end


% [centra_p1,ind1]=sort(centra(:,1));
% xp1=[];
% yprobp1=[];
% for i=1:size(centra_p1)
%     xp1=[xp1 centra_p1(i)-deltap(ind1(i),1) centra_p1(i)-deltap(ind1(i),1)];
%     yprobp1=[yprobp1 prob(ind1(i)) prob(ind1(i))];
% end
% plot((1:41),datay(:,1),(1:41),obsy)
% plot(datax(:,1),datay(:,1),ToutX,OutX(:,1))
% 
% N=1000;
% filename = 'theta_sample.txt';
% [theta1,delimiterOut]=importdata(filename);
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


%project 2d to 1d
ncell=size(prob2,1);
den=zeros(ncell,2);
for i=1:ncell
    den(i,1)=exp(prob2(i,1))/deltap(i,1)/(2^1);%lower bound
    den(i,2)=exp(prob2(i,2))/deltap(i,1)/(2^1);%upper bound
end

centra_1=[centra(1,1)];
den_1= [den(1,1) den(1,2)];
deltap_1=[deltap(1,1)];
for i=2:ncell
    ind_temp1=find(centra_1(:,1)==centra(i,1));
    if size(ind_temp1,1)==0
        
        centra_1=[centra_1;centra(i,1)];
        den_1=[den_1; den(i,1) den(i,2)];
        deltap_1=[deltap_1; deltap(i,1)];
    else           
        ind_den_1=ind_temp1;
        den_1(ind_den_1,:)=den_1(ind_den_1,:)+[den(i,1) den(i,2)];
       
    end

end

deInd=[];
for i=1:size(centra_1(:,1),1)
    for j=1:size(centra_1(:,1),1)
        if (i~=j)
            if(centra_1(i,1)+deltap_1(i,1) <=centra_1(j,1)+deltap_1(j,1) && centra_1(i,1)-deltap_1(i,1) >= centra_1(j,1)-deltap_1(j,1))
                den_1(i,:)=den_1(i,:)+den_1(j,:);
                deInd=[deInd;j];
            end
            
        end
    end
end


%plot(centra_1(:,1),den_1(:,1),'.b',centra_1(:,1),den_1(:,2),'.r')
deInd_uniq = unique(deInd);
centra_1new=centra_1(setdiff(1:end,deInd_uniq),1);
den_1lo=den_1(setdiff(1:end,deInd_uniq),1);
den_1up=den_1(setdiff(1:end,deInd_uniq),2);
den_1new=[den_1lo den_1up];
subplot(2,1,1);
%plot(centra_1new(:,1),den_1new(:,1),'.b',centra_1new(:,1),den_1new(:,2),'.r')
nCnew=size(centra_1new,1);
Sden_1new=zeros(nCnew,2);
Scentra_1new=zeros(nCnew,1);
[sortV , sortIndex] = sort(centra_1new(:,1));
Scentra_1new(:,1)=sortV;
for i=1:nCnew
    Sden_1new(i,1) = den_1new( sortIndex(i), 1);
end
stairs(Scentra_1new(:,1),Sden_1new(:,1),'*-b','LineWidth',1) 
hold on;
for i=1:nCnew
    Sden_1new(i,2) = den_1new( sortIndex(i), 2);
end
stairs(Scentra_1new(:,1),Sden_1new(:,2),'*-r','LineWidth',1)


ncell=size(prob2,1);
den=zeros(ncell,2);
for i=1:ncell
    den(i,1)=exp(prob2(i,1))/deltap(i,2)/(2^1);%lower bound
    den(i,2)=exp(prob2(i,2))/deltap(i,2)/(2^1);%upper bound
end

centra_2=[centra(1,2)];
den_2= [den(1,1) den(1,2)];
deltap_2=[deltap(1,2)];
for i=2:ncell
    ind_temp1=find(centra_2(:,1)==centra(i,2));
    if size(ind_temp1,1)==0
        
        centra_2=[centra_2;centra(i,2)];
        den_2=[den_2; den(i,1) den(i,2)];
        deltap_2=[deltap_2; deltap(i,2)];
    else           
        ind_den_2=ind_temp1;
        den_2(ind_den_2,:)=den_2(ind_den_2,:)+[den(i,1) den(i,2)];
       
    end

end

deInd=[];
for i=1:size(centra_2(:,1),1)
    for j=1:size(centra_2(:,1),1)
        if (i~=j)
            if(centra_2(i,1)+deltap_2(i,1) <=centra_2(j,1)+deltap_2(j,1) && centra_2(i,1)-deltap_2(i,1) >= centra_2(j,1)-deltap_2(j,1))
                den_2(i,:)=den_2(i,:)+den_2(j,:);
                deInd=[deInd;j];
            end
            
        end
    end
end

deInd_uniq = unique(deInd);
%setdiff(1:end,deInd_uniq)
centra_2new=centra_2(setdiff(1:end,deInd_uniq),1);
den_2lo=den_2(setdiff(1:end,deInd_uniq),1);
den_2up=den_2(setdiff(1:end,deInd_uniq),2);
den_2new=[den_2lo den_2up];


subplot(2,1,2);
%plot(centra_2new(:,1),den_2new(:,1),'.b',centra_2new(:,1),den_2new(:,2),'.r')
nCnew=size(centra_2new,1);
Sden_2new=zeros(nCnew,2);
Scentra_2new=zeros(nCnew,1);
[sortV , sortIndex] = sort(centra_2new(:,1));
Scentra_2new(:,1)=sortV;
for i=1:nCnew
    Sden_2new(i,1) = den_2new( sortIndex(i), 1);
end
stairs(Scentra_2new(:,1),Sden_2new(:,1),'*-b','LineWidth',1) 
hold on;
for i=1:nCnew
    Sden_2new(i,2) = den_2new( sortIndex(i), 2);
end
stairs(Scentra_2new(:,1),Sden_2new(:,2),'*-r','LineWidth',1)
