function [YhatDL, meanY, GammaY, likeGammaY, likethetaY] = estim(Y, gammaY)

T = size(Y,1);
meanY = mean(Y);

rhoY = gammaY/gammaY(1); 

% disp('for GammaY')
% GammaY = zeros(T,T);
% for h = 0:T-1
%     GammaTemp = gammaY(h+1)*ones(T-h,1);
%     GammaY = GammaY + diag(GammaTemp,h) + diag(GammaTemp,-h);
% end
% 
% disp('for loglikeGammaY')
% %Log likehood for GammaY
% likeGammaY = ((2*pi())^(-T/2))*(det(GammaY)^(-1/2))*exp((-1/2)*(Y'*inv(GammaY)*Y));

% %Durbin Levinson
% disp('for DL')
% aDL = zeros(T,T);
% vDL = zeros(T,1);
% YhatDL = zeros(T,1);
% aDL(1,1) = gammaY(2)/gammaY(1);
% vDL(1,1) = gammaY(1);
% YhatDL(1,1) = Y(1);
% for t = 2:T
%     DLTemp = 0;
%     for j = 1:t-1
%         DLTemp = DLTemp + aDL(j,t-1)*gammaY(t-j);
%     end
%     aDL(t,t) = (gammaY(t)-DLTemp)/vDL(t-1,1);
%     aDL(1:t-1,t) = aDL(1:t-1,t-1) - aDL(t,t)*aDL(t-1:-1:1,t-1);
%     vDL(t,1) = vDL(t-1,1)*(1-aDL(t,t)^2);
%     YhatDL(t,1) = aDL(:,t)'*Y;
% end

%Innovations
disp('for Inn')
thetaInn = zeros(T,T);
vInn = zeros(T,1);
YhatInn = zeros(T,1);
vInn(1,1) = gammaY(1);
thetaInn(1,1) = gammaY(2)/gammaY(1);
for t = 2:T
    InnTemp = 0;
    for k = 2:t-1
            for j = k-1:-1:0
            InnTemp = InnTemp + thetaInn(k-j,k)*thetaInn(t-j,t)*vInn(j+1,1);
            end
        thetaInn(t-k,t) = (gammaY(t-k)-InnTemp)/vInn(t-k);
    end
    InnTemp = 0;
    Yhattemp = 0;
    for l = 2:t-1
        InnTemp = InnTemp + (thetaInn(t-l,t)^2)*vInn(l-1);
        Yhattemp = Yhattemp + thetaInn(l-1,t)*(Y(t+1-l) - YhatInn(t+1-j));
    end
    vInn(t,1) = gammaY(1) - InnTemp;
    YhatInn(t+1) = Yhattemp;
    disp(t)
end

disp('end');