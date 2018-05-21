clear all
close all
clc
disp('MA1')
rng('default')     

% Monte-Carlo runs 
runs = 10;

%AR(1) model specs
n = 300;
Y   = NaN(n,runs);
theta = [0 ; 0.2 ; -0.5]; %[mu_y sigma phi_1]
p = 1;
q = 0;

%residuals
epsY = theta(2)*randn(n,runs);

%Generate the AR(1) process
Y(1,:) = epsY(1,:);
for t = 1:n-1
   Y(t+1,:) = theta(1) + epsY(t+1,:) + theta(3)*epsY(t,:);
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
thetaStart = [0.1 ; 0.1 ; -0.7]; 
options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

%% DL/Inn and MLE

for i = 1:runs
    
     objfun = @(thetaStart)(-loglikeMA1(Y(:,i), thetaStart, n));
     [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
%      disp(i)

% ================

    % finding all sample autocov values
%     meanY = mean(Y(:,i));
%     gammaY = zeros(n,1);
%     for k = 1:n
%         for t = 1:n+1-k
%             gammaY(k) = gammaY(k) + (Y(t,i)-meanY)*(Y(t+k-1,i)-meanY);
%         end
%         gammaY(k) = gammaY(k)/(n+1-k);
%     end 
    
%     [YhatDL(:,i), vDL(:,i)] = durblev(Y(:,i), gammaY);
%     [YhatInn(:,i), vInn(:.i)] = innov(Y, gammaY);

% % DL
%         objfun = @(thetaStart)(durblev(Y(:,i)-thetaStart(1),thetaStart));
% %Inn        
% %         objfun = @(thetaStart)(innov(Y(:,i),gammaY,thetaStart));
% %MLE    
%         [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);


    if mod(i,100)==0
        disp(i);
    end
end
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
% display(mean(theta_mle(1,:)));
display(mean(theta_mle(2,:)));
display(mean(theta_mle(3,:)));