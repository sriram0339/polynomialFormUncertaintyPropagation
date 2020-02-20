% Simulate the Rimless wheel model
nSims = 100000;
nSteps = 5000;
xVals = zeros(nSims,1);
for i = 1:nSims
    x = 0.2 * rand - 0.1;
    theta = pi/6;
    g = 10.0;
    for j = 1: nSteps
       w = pi/180 *(8 + 1.5 * randn);
       beta1 = theta/2 + w;
       beta2 = theta/2 - w;
       x = cos(theta)^2 * ( x + 2*g*(1- cos(beta1))) - 2*g * (1-cos(beta2));
    end
    xVals(i,1) = x;
end

figure(1);
hist(xVals,100);
fprintf(1, 'x: expectation = %f, range = [%f, %f] \n', mean(xVals), min(xVals), max(xVals));