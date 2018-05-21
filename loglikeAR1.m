function [loglike] = loglikeAR1(Y, thetaStart, T)

Yhat = zeros(T,1);
Yhat(2:T) = thetaStart(3)*Y(1:T-1);
v = thetaStart(2)^2;

objfunTemp = 0;
    for t = 2:T
        objfunTemp = objfunTemp + log(v) + (Y(t) - Yhat(t))^2/v;
    end
loglike = -T/(2*pi()) - 1/2*objfunTemp;

end