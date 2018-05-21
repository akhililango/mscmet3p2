function [loglike] = loglikeARMA11(Y, thetaStart, T)

epshat = zeros(T,1);
epshat(1) = Y(1);
for t =2:T
    epshat(t) = Y(t) - thetaStart(3)*Y(t-1) - thetaStart(4)*epshat(t-1);
end
v = thetaStart(2)^2;




objfunTemp = 0;
    for t = 2:T
        objfunTemp = objfunTemp + log(v) + (epshat(t))^2/v;
    end
loglike = -T/(2*pi()) - 1/2*objfunTemp;

end