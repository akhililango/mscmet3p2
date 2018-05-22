clear all
close all
clc
disp('AR1')
rng('default')     

% global Y epsY
% Monte-Carlo runs 
runs = 100;

%AR(1) model specs
T = 300;
Y   = NaN(T,runs);
theta = [0.2 ; 0.7]; %[mu_y sigma phi_1]
p = 1;
q = 0;

%residuals
epsY = theta(1)*randn(T,runs);

%Generate the AR(1) process
Y(1,:) = epsY(1,:);
for t = 1:T-1
   Y(t+1,:) = theta(2)*Y(t,:) + epsY(t+1,:);  
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
thetaStart = [0.1 ; 0.5]; 
options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

%% DL/Inn and MLE

for i = 1:runs
    
     objfun = @(thetaStart)(-loglikeAR1(Y, thetaStart, T));
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
%          objfun = @(thetaStart)(durblev(Y(:,i)-thetaStart(1),thetaStart));
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

% fprintf ('Log Likelihood value = %g \r', durblev());
theta_mle_1 = mean(theta_mle(1,:));
theta_mle_2 = mean(theta_mle(2,:));
display(theta_mle_1);
display(theta_mle_2);

SE_21 = theta_mle_2 + 1.96*theta_mle_1/sqrt(T);
SE_22 = theta_mle_2 - 1.96*theta_mle_1/sqrt(T);

f1 = figure;
histfit(theta_mle(2,:),25,'kernel');
line([theta_mle_2, theta_mle_2], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([theta_mle_2+1.96*theta_mle_1/sqrt(T) theta_mle_2+1.96*theta_mle_1/sqrt(T) NaN theta_mle_2-1.96*theta_mle_1/sqrt(T) theta_mle_2-1.96*theta_mle_1/sqrt(T)] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')

%% Criteria Testing

    Y_AIC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 1);
    Y_BIC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 2);
    Y_AICC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 3);
    
    Y_AIC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 1);
    Y_BIC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 2);
    Y_AICC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 3);
