function [loglike] = loglikeMA1(Y, thetaStart, T)

epshat = zeros(T,1);
epshat(1) = Y(1) - mean(Y);
epshat(2:T) = Y(2:T) - thetaStart(3)*epshat(1:T-1);
v = thetaStart(2)^2;




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
        objfunTemp = objfunTemp + log(v) + (Y(t) - thetaStart(3)*epshat(t-1))^2/v;
 
%         objfunTemp = objfunTemp + log(v) + (Y(t) - Yhat(t))^2/v;
    end
    
    
    
loglike = -T/(2*pi()) - 1/2*objfunTemp;

end