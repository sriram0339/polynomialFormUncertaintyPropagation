nSims = 100;
nSteps = 20;
xVals = zeros(nSims,nSteps);
yVals = zeros(nSims,nSteps);
vVals = zeros(nSims,nSteps);
psiVals = zeros(nSims,nSteps);
for i = 1:nSims
    x = rand * 0.2 - 0.1;
    y = rand * 0.2 - 0.5;
    v = rand * 1.5 + 6.5;
    psi = 0.1 + 0.1 * randn;
    dpsi = randn * 0.1;
    for j = 1:nSteps
        
        xVals(i,j) = x;
        yVals(i,j) = y;
        vVals(i,j) = v;
        psiVals(i,j )= psi;
        x = x + 0.1 * v * cos(psi);
        y = y + 0.1 * v * sin(psi);
        v = v + 0.1 * (-0.5 *  (v - 10) + 0.2 * rand - 0.1);
        psi =  0.1 * randn + 0.1;
           
    end
   figure(1);
   hold on;
   plot(xVals(i,:), yVals(i,:),'b');
   

end
 figure(2);
 hold on;
 scatter(xVals(:,10), yVals(:,10),'.','k');
 scatter(xVals(:,5), yVals(:,5),'.','k');
 scatter(xVals(:,1), yVals(:,1),'.','k');
 scatter(xVals(:,15), yVals(:,15),'.','k');
 scatter(xVals(:,20), yVals(:,20),'.','k');
 fprintf(1, 'x(10): range [%f, %f], expectation : %f ', ...
     min(xVals(:,11)), max(xVals(:,11)), mean(xVals(:,11)));
  fprintf(1, 'y(10): range [%f, %f], expectation : %f ', ...
     min(yVals(:,11)), max(yVals(:,11)), mean(yVals(:,11)));
