d=2;
maxtheta=[0.4; 0.2];
mintheta=[0.2; 0.1];
minw=-(maxtheta-mintheta)/2;
maxw=(maxtheta-mintheta)/2;


is_grad=1;
is_symbolic=1;
% number of intervals on each dimension in variational histogram
n_int=20;
syms w2 w3
w = sym([w2 w3]);