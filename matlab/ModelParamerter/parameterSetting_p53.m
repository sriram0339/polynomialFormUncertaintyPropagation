d=2;
maxtheta=[0.51; 0.51];
mintheta=[0.49; 0.49];
minw=-(maxtheta-mintheta)/2;%-10;
maxw=(maxtheta-mintheta)/2;%10;

is_grad=0;%1;
is_symbolic=1;
% number of intervals on each dimension in variational histogram
n_int=50;%20;
%n_int_MG=10000;
syms w6 w7
w = sym([w6 w7]);