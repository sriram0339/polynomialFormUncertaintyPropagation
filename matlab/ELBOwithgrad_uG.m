%function [f,g] = ELBOwithgrad_uG(theta,expectation_q,negativeH,gradF2,gradNegativeH)
function [f,g] = ELBOwithgrad_uG(theta,sym_theta,f,negativeH,gradF2g,gradNegativeH)
%global expectation_q negativeH gradF2 gradNegativeH
% Calculate objective f
f = -double(vpa(subs(f,sym_theta,theta)))+ negativeH(theta);%;-expectation_q(theta) + negativeH(theta);
if nargout > 1 % gradient required
    g =-double(vpa(subs(gradF2g,sym_theta,theta))) + gradNegativeH(theta);% -gradF2(theta) + gradNegativeH(theta);
end

