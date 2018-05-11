function [loglikethetaY] = loglikeARMA(thetaStart, Y, epsY, p, q)
iT = size(Y,1);

sigma = thetaStart(p+q+2);
Yhat = 0;

Ytilde(1:p) = Y(1);

    for t = p:iT-1
            Yhat = 0;
            if p ~= 0
                for j = 1:p 
                    Yhat = Yhat + thetaStart(j+1)*Y(t+1-j);
                end
            end
            if q~=0
                for j = p+1:p+q 
                    Yhat = Yhat + thetaStart(j+1)*epsY(t+p+1-j);
                end
            end
            Ytilde(t+1) = Y(t+1) - Yhat;
    end
    
loglikethetaY = log(((2*pi()*sigma^2)^(-iT/2))*exp((-1/(2*sigma^2))*sum((Ytilde).^2)));

end
