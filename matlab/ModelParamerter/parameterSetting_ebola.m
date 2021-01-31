d=3;
maxtheta=[0.5; 0.35; 0.35];
mintheta=[0.25; 0.25; 0.25];

minw=-(maxtheta-mintheta)/2;%-10;
maxw=(maxtheta-mintheta)/2;%10;


is_grad=1;
is_symbolic=1;
% number of intervals on each dimension in variational histogram
n_int=40;%20;
syms w5 w6 w7
w = sym([w5 w6 w7]);

% for the general factorization
size_C=2;
Clusters{1}=[1,3];
Clusters{2}=[2];