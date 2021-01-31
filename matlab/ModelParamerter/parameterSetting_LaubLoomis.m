d=13;
maxtheta=[3.2; 3; 2.5; 5; 3.5; 5; 3.8; 3.8; 2.9; 2.9; 3.3; 6.4; 4.7];
mintheta=[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
minw=-(maxtheta-mintheta)/2;%-10;
maxw=(maxtheta-mintheta)/2;%10;


is_grad=0;
is_symbolic=0;
% number of intervals on each dimension in variational histogram
n_int=5;%10;
syms w7 w8 w9 w10 w11 w12 w13 w14 w15 w16 w17 w18 w19
w = sym([w7 w8 w9 w10 w11 w12 w13 w14 w15 w16 w17 w18 w19]);

% for the general factorization
size_C=5;
Clusters{1}=[10];
Clusters{2}=[5 6 11];
Clusters{3}=[1 2 7];
Clusters{4} =[4 8 12];
Clusters{5} = [3 9 13];