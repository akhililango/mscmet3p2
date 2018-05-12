function [YhatInn, vInn, meanY] = innov(Y, gammaY)

T = size(Y,1);
meanY = mean(Y);

%Innovations
disp('for Inn')
thetaInn = zeros(T,T);
vInn = zeros(T,1);
YhatInn = zeros(T,1);
vInn(1,1) = gammaY(1);
thetaInn(1,1) = 0.6180;

for t = 2:T
    Yhattemp = 0;
    for j = 1:t-1
        Yhattemp = Yhattemp + thetaInn(j,t-1)*(Y(t-j) - YhatInn(t-j));
    end
    YhatInn(t) = Yhattemp;
    
    thetaInnTemp = 0;
%     vInnTemp = 0;
    thetaInn(t,t) = (gammaY(t)-thetaInn(t-1,t-1)^2*vInn(t-1,1))/vInn(t-1,1);
    for k = 1:t-1
            for j = 0:k-1
            thetaInnTemp = thetaInnTemp + thetaInn(k+1-j,k+1)*thetaInn(t-j,t)*vInn(j+1,1);
%             vInnTemp = vInnTemp + (thetaInn(t-j,t)^2)*vInn(j+1);
            end
        thetaInn(t-k,t) = (gammaY(t-k)-thetaInnTemp)/vInn(k+1,1);
%         vInn(t,1) = gammaY(1) - vInnTemp;
    end
    
    vInnTemp = 0;
    for j = 0:t-2
        vInnTemp = vInnTemp + (thetaInn(t-j,t-1)^2)*vInn(j+1);
    end
    vInn(t,1) = gammaY(1) - vInnTemp;
    disp(t)
end

disp('end');