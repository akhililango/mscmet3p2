function [gammaY] = gammaAR1(T, thetaStart)

psi = thetaStart(3);
%to find gammaY
gammaY = zeros(T,1);
for h = 1:T
    gammaY(h) = psi^h/(1-psi^2);
end
disp('');