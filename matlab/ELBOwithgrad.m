function [f,g] = ELBOwithgrad(theta,expectation_q,negativeH,gradF2,gradNegativeH)
%function [f,g] = ELBOwithgrad(theta)
%global expectation_q negativeH gradF2 gradNegativeH
% Calculate objective f
f = -expectation_q(theta) + negativeH(theta);
if nargout > 1 % gradient required
    g = -gradF2(theta) + gradNegativeH(theta);
end

