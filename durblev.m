function [loglike, YhatDL, vDL] = durblev(Y, thetaStart)

T = size(Y,1);
rng('default')

gammaY = gammaMA1(T, thetaStart);

%Durbin Levinson
% disp('for DL')
aDL = zeros(T,T);
vDL = zeros(T,1);
YhatDL = zeros(T,1);
% aDL(1,1) = thetaStart(1);
aDL(2,2) = gammaY(2)/gammaY(1);
vDL(1,1) = gammaY(1);
vDL(2,1) = gammaY(1)-gammaY(2);
YhatDL(1,1) = aDL(1,1);
YhatDL(2,1) = aDL(2,2)*Y(1);
for t = 3:T
    DLTemp = aDL(:,t-1)'*ones(T,1);
    aDL(t,t) = (gammaY(t)-DLTemp)/vDL(t-1,1);
%     aDL(t,t) = min(aDL(t,t),0.99);
%     aDL(t,t) = max(aDL(t,t),-0.99);
    aDL(2:t-1,t) = aDL(2:t-1,t-1) - aDL(t,t)*aDL(t-1:-1:2,t-1);
    vDL(t,1) = vDL(t-1,1)*(1-aDL(t,t)^2);
    YhatDL(t,1) = aDL(2:t,t)'*Y(1:t-1);
end

objfunTemp = 0;
    for t = 1:T
        if vDL(t)>0
        objfunTemp = objfunTemp + log(vDL(t)) + (Y(t) - YhatDL(t))^2/vDL(t);
        end
    end
loglike = -(-T/(2*pi()) - 1/2*objfunTemp);
