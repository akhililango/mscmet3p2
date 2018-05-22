function [loglike] = loglikeAR1(Y, thetaStart, T, criteria)

if nargin < 4
    criteria = 0;
end

Yhat = zeros(T,1);
Yhat(2:T) = thetaStart(2)*Y(1:T-1);
v = thetaStart(1)^2;

objfunTemp = 0;
    for t = 2:T
        objfunTemp = objfunTemp + log(v) + (Y(t) - Yhat(t))^2/v;
    end
    
    if criteria == 0
        loglike = -T/(2*pi()) - 1/2*objfunTemp;
    elseif criteria == 1
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 4;         %AIC
    elseif criteria == 2
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 4*T/(T-3);         %BIC
    elseif criteria == 3
        loglike = -2*(-T/(2*pi()) - 1/2*objfunTemp) + 2*log(T)/T;         %AICC
    else
        disp('Invalid criteria entry')
    end

end