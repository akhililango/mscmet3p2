function [gammaY] = gammaMA1(T, thetaStart)

%to find gammaY
gammaY = zeros(T,1);
gammaY(1) = (1+thetaStart(3)^2)*thetaStart(2)^2;
gammaY(2) = thetaStart(3)*thetaStart(2)^2;
disp('');