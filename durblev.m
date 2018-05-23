function [like] = durblev(Yreal, thetaStart)

T = size(Yreal,1);
% rng('default')

%%
% %residuals
% epsY = thetaStart(1)*randn(T,1);
% 
% Y = zeros(T,1);
% %Generate the AR(1) process
% Y(1) = Yreal(1);
% for t = 1:T-1
%    Y(t+1) = thetaStart(1) + thetaStart(2)*Y(t) + epsY(t+1);  
% end
% 
% rhoY = autocorr(Y,T-1);
% gammaY = abs(rhoY*var(Y));

%%
gammaY = gammaAR1(T, thetaStart);

%% Durbin Levinson
% disp('for DL')
aDL = zeros(T,T);
vDL = zeros(T,1);
YhatDL = zeros(T,1);
% aDL(1,1) = thetaStart(1);
aDL(2,2) = gammaY(2)/gammaY(1);
vDL(1,1) = gammaY(1);
vDL(2,1) = gammaY(1)*(1-gammaY(2)^2/gammaY(1)^2);
YhatDL(1,1) = Yreal(1);
YhatDL(2,1) = aDL(2,2)*Yreal(1);
for t = 3:T
%     DLTemp = aDL(:,t-1)'*ones(T,1);
    DLTemp = 0;
    for j = 1:t-1
        DLTemp = DLTemp + aDL(j,t-1)*gammaY(t-j);
    end
    aDL(t,t) = (gammaY(t)-DLTemp)/vDL(t-1,1);
%     aDL(t,t) = min(aDL(t,t),0.99);
%     aDL(t,t) = max(aDL(t,t),-0.99);
    aDL(2:t-1,t) = aDL(2:t-1,t-1) - aDL(t,t)*aDL(t-1:-1:2,t-1);
    vDL(t,1) = vDL(t-1,1)*(1-aDL(t,t)^2);
    YhatDL(t,1) = aDL(2:t,t)'*Yreal(1:t-1,1);
end


%% Likelihood
objfunTemp = 0;
    for t = 1:T
        objfunTemp = objfunTemp + (Yreal(t) - YhatDL(t))^2/vDL(t);
    end
% loglike = -(-T/(2*pi()) - 1/2*objfunTemp);

like =  (-0.5*sum(log(vDL))-objfunTemp/2); 

% M = toeplitz(vDL);
% S = Yreal'*M*Yreal;
% g = -log(det(M));
% like = -(-(T/2)*log(S/T)-0.5*g);
disp('');