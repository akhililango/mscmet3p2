function [loglike, YhatInn, vInn] = innov(Y, gammaY, thetaStart)

n = size(Y,1);
meanY = mean(Y);

%Innovations
disp('for Inn')
thetaInn = zeros(n,n);
vInn = zeros(n,1);
YhatInn = zeros(n,1);
vInn(1,1) = gammaY(1);
thetaInn(1,1) = thetaStart(2);

for t = 2:n
    Yhattemp = 0;
    for j = 1:t-1
        Yhattemp = Yhattemp + thetaInn(j,t-1)*(Y(t-j) - YhatInn(t-j));
    end
    YhatInn(t) = Yhattemp; 
    
    thetaInnTemp = 0;
%     vInnTemp = 0;
    for k = 1:t-1
            for j = 0:k-1
            thetaInnTemp = thetaInnTemp + thetaInn(k-j,k)*thetaInn(t-1-j,t-1)*vInn(j+1,1);
%             vInnTemp = vInnTemp + (thetaInn(t-j,t)^2)*vInn(j+1);
            end
        thetaInn(t-k,t) = (gammaY(t-k)-thetaInnTemp)/vInn(k,1);
%         vInn(t,1) = gammaY(1) - vInnTemp;
    end
    
    vInnTemp = 0;
    for j = 0:t-2
        vInnTemp = vInnTemp + (thetaInn(t-1-j,t-1)^2)*vInn(j+1);
    end
    vInn(t,1) = gammaY(1) - vInnTemp;
    disp(t)
end

objfunTemp = 0;
    for t = 1:n
        objfunTemp = objfunTemp + log(vInn(t)) + (Y(t) - YhatInn(t))^2/vInn(t);
    end
loglike = -n/2/pi() - 1/2*objfunTemp;

disp('end');
