d=5;
maxtheta=[0.0011; 0.0011; 0.5; 1.0; 1.0];
mintheta=[0.0009; 0.0009; 0.0; 0.0; 0.0];
minw=-(maxtheta-mintheta)/2;%-10;
maxw=(maxtheta-mintheta)/2;%10;


is_grad=0;
is_symbolic=1;
% number of intervals on each dimension in variational histogram
n_int=10;
syms w5 w6 w7 w8 w9
w = sym([w5 w6 w7 w8 w9]);