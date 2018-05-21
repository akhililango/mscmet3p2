clear all
close all
clc
disp('AR1')
rng('default')     

% Monte-Carlo runs 
runs = 1000;

%AR(1) model specs
n = 300;
Y   = NaN(n,runs);
theta = [0 ; 0.9 ; 0.2]; %[mu_y phi_1 sigma]
p = 1;
q = 0;

%residuals
epsY = theta(3)*randn(n,runs);

%Generate the AR(1) process
Y(1,:) = epsY(1,:);
for t = 1:n-1
   Y(t+1,:) = theta(1) + theta(2)*Y(t,:) + epsY(t+1,:);
   
end

% % To plot some stuff 
% subplot(2,1,1);   
% autocorr(Y(:,1),100);
% 
% subplot(2,1,2);
% [f,xi] = ksdensity(Y(:,1));
% ts1 = f;
% plot(ts1);


% estimation in ARMA(p,q) model 
thetaStart = [0.3 ; 0.5 ; 0.15]; 
options = optimset('TolX', 0.0001, 'Display', 'off', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

%% DL/Inn and MLE

for i = 1:runs
    
%      objfun = @(thetaStart)(-loglikeARMA(thetaStart, Y(:,i), epsY(:,i), p, q));
%      [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
%      disp(i)

% ================

    % finding all sample autocov values
    meanY = mean(Y(:,i));
    gammaY = zeros(n,1);
    for k = 1:n
        for t = 1:n+1-k
            gammaY(k) = gammaY(k) + (Y(t,i)-meanY)*(Y(t+k-1,i)-meanY);
        end
        gammaY(k) = gammaY(k)/(n+1-k);
    end 
    
%     [YhatDL(:,i), vDL(:,i)] = durblev(Y(:,i), gammaY);
%     [YhatInn(:,i), vInn(:.i)] = innov(Y, gammaY);

% DL
%         objfun = @(thetaStart)(durblev(Y(:,i),gammaY,thetaStart));
%Inn        
        objfun = @(thetaStart)(innov(Y(:,i),gammaY,thetaStart));
%MLE    
        [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
    disp(i)

end

vDL = abs(vDL);

%% MLE

% for i = 1:runs
% %     objfun = @(thetaStart)(-loglikeARMA(thetaStart, Y(:,i), epsY, p, q));
% %     objfunTemp = 0;
% %     for t = 1:n
% %         objfunTemp = objfunTemp + log(vDL(t,i)) + (Y(t,i) - YhatDL(t,i))^2/vDL(t,i);
% %     end
% %     objfun = @(thetaStart)(durblev(Y));
%         
%     [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
%     disp(i)
% end

%% Display

fprintf ('Log Likelihood value = %g \r', -dLogLik*runs);
display(mean(theta_mle(1,:)));
display(mean(theta_mle(2,:)));
display(mean(theta_mle(3,:)));
histfit(theta_mle(2,:));
