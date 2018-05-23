function [loglike] = loglikeMA1(Y, thetaStart, T, criteria)

if nargin < 4
    criteria = 0;
end

epshat = zeros(T,1);
epshat(1) = Y(1);% - mean(Y);
epshat(2:T) = Y(2:T) - thetaStart(2)*epshat(1:T-1);
v = thetaStart(1)^2;

% 
% Yhat = zeros(T,1);
% for t = 2:T
%     Yhat(t) = mean(Y);
%     for j = 1: t
%         Yhat(t) = Yhat(t) + (-1)^(j+1)*(thetaStart(3))^(j)*(Y(t+1-j) - mean(Y)); 
%     end
% end
% v = thetaStart(2)^2;



objfunTemp = 0;
    for t = 2:T
        objfunTemp = objfunTemp + log(v) + (Y(t) - thetaStart(2)*epshat(t-1))^2/v;
 
%         objfunTemp = objfunTemp + log(v) + (Y(t) - Yhat(t))^2/v;
    end
    
% loglike = -T/(2*pi()) - 1/2*objfunTemp;

    if criteria == 0
%         loglike = -T/(2*pi()) - 1/2*objfunTemp;
        loglike = - 1/2*objfunTemp;
    elseif criteria == 1
        loglike = 2*(-T/(2*pi()) - 1/2*objfunTemp) + 4;         %AIC
    elseif criteria == 2
        loglike = 2*(-T/(2*pi()) - 1/2*objfunTemp) + 4*T/(T-3);         %BIC
    elseif criteria == 3
        loglike = 2*(-T/(2*pi()) - 1/2*objfunTemp) + 2*log(T)/T;         %AICC
    else
        disp('Invalid criteria entry')
    end

end