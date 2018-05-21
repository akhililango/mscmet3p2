function [loglike, YhatDL, vDL] = durblev(Y, thetaStart, p, q)

T = size(Y,1);
rng('default')
%residuals
epsY = thetaStart(2)*randn(T);

%Generate the AR(1) process
y(1) = Y(1);
for t = 1:T-1
   y(t+1) = thetaStart(1) + thetaStart(3)*y(t) + epsY(t+1);  
end

rhoY = autocorr(y,T-1);
gammaY = abs(rhoY*var(y));

% disp('for GammaY')
% GammaY = zeros(T,T);
% for h = 0:T-1
%     GammaTemp = gammaY(h+1)*ones(T-h,1);
%     GammaY = GammaY + diag(GammaTemp,h) + diag(GammaTemp,-h);
% end

% disp('for loglikeGammaY')
% %Log likehood for GammaY
% likeGammaY = ((2*pi())^(-T/2))*(det(GammaY)^(-1/2))*exp((-1/2)*(Y'*inv(GammaY)*Y));

%Durbin Levinson
disp('for DL')
aDL = zeros(T,T);
vDL = zeros(T,1);
YhatDL = zeros(T,1);
aDL(1,1) = gammaY(2)/gammaY(1);
vDL(1,1) = gammaY(1);
YhatDL(1,1) = y(1);
for t = 2:T
    DLTemp = 0;
    for j = 1:t-1
        DLTemp = DLTemp + aDL(j,t-1)*gammaY(t-j);
    end
    aDL(t,t) = (gammaY(t)-DLTemp)/vDL(t-1,1);
    aDL(1:t-1,t) = aDL(1:t-1,t-1) - aDL(t,t)*aDL(t-1:-1:1,t-1);
    vDL(t,1) = vDL(t-1,1)*(1-aDL(t,t)^2);
    YhatDL(t,1) = aDL(:,t)'*Y;
end

objfunTemp = 0;
    for t = 1:T-1
        objfunTemp = objfunTemp + log(vDL(t)) + (y(t+1) - YhatDL(t+1))^2/vDL(t);
    end
loglike = -T/2/pi() - 1/2*objfunTemp;
    
disp('end')
