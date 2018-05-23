function [loglike] = loglikeARMA11(Y, thetaStart, T, criteria)

if nargin < 4
    criteria = 0;
end

epshat = zeros(T,1);
epshat(1) = Y(1);
for t =2:T
    epshat(t) = Y(t) - thetaStart(2)*Y(t-1) - thetaStart(3)*epshat(t-1);
end
v = thetaStart(1)^2;

objfunTemp = 0;
    for t = 2:T
        objfunTemp = objfunTemp + log(v) + (Y(t) - thetaStart(2)*Y(t-1) - thetaStart(3)*epshat(t-1))^2/v;
    end
    
loglike = -T/(2*pi()) - 1/2*objfunTemp;

    if criteria == 0
        loglike = -T/(2*pi()) - 1/2*objfunTemp;
    elseif criteria == 1
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 6;         %AIC
    elseif criteria == 2
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 6*T/(T-3);         %BIC
    elseif criteria == 3
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 3*log(T)/T;         %AICC
    else
        disp('Invalid criteria entry')
    end

end