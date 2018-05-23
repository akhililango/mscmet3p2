% function [] = modelestimation(process)

clear all
close all
clc
disp('AR1')
rng('default')     

% global Y epsY
% Monte-Carlo runs 
runs = 100;

%ARMA(1,1) model specs
T = 300;
Y   = NaN(T,runs);
theta = [0.4 ; 0.3; -0.4]; 
p = 1;
q = 1;

%residuals
epsY = theta(1)*randn(T,runs);

%Generate the ARMA(1,1) process
Y(1,:) = epsY(1);
for t = 1:T-1
   Y(t+1,:) = theta(2)*Y(t,:) + theta(3)*epsY(t,:) + epsY(t+1,:);
   
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

for i = 1:80
    
     objfun = @(thetaStart)(-loglikeAR1(Y(:,i), thetaStart, T));
     [theta_mle_AR1(:,i), dLogLikAR1(i),~,~,~,hess] = fminunc(objfun, thetaStart, options);
%      disp(i)
     invhess = inv(hess);
     SEAR1(i) = 1.96*sqrt(invhess(2,2));
    
     objfun = @(thetaStart)(-loglikeMA1(Y(:,i), thetaStart, T));
     [theta_mle_MA1(:,i), dLogLikMA1(i),~,~,~,hess] = fminunc(objfun, thetaStart, options);
     invhess = inv(hess);
     SEMA1(i) = 1.96*invhess(2,2);
% 
%      objfun = @(thetaStart)(-loglikeARMA11(Y(:,i), [0.1 ; 0.5 ; 0.1], T));
%      [theta_mle_ARMA11(:,i), dLogLikARMA11(i),~,~,~,hess] = fminunc(objfun, [0.1 ; 0.5 ; 0.1], options);
%      invhess = inv(hess);
%      SEARMA11 = 1.96*invhess(2,2);
%      SEARMA11th = 1.96*invhess(3,3);
% ================
    
    res(i) = Hevia_arma_mle(Y(:,i), 1, 1);
    theta_mle_ARMA11(1,i) = res(i).sigma;
    theta_mle_ARMA11(2,i) = 2*normcdf(res(i).ar)-1;
    theta_mle_ARMA11(3,i) = 2*normcdf(res(i).ma)-1;
    invhess = inv(res(i).hess);
    SEARMA11 = 1.96*invhess(1,1);
    SEARMA11th = 1.96*invhess(2,2);
    dLogLikARMA11(i) = res(i).loglike;

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
theta_mle_AR1_1 = mean(theta_mle_AR1(1,:));
% theta_mle_AR1_2 = mean((2*normcdf(theta_mle_AR1(2,:))-1));
theta_mle_AR1_2 = mean(theta_mle_AR1(2,:));
loglikeAR1 = mean(dLogLikAR1);
display(theta_mle_AR1_1);
display(theta_mle_AR1_2);

% SE_21 = theta_mle_2 + 1.96*theta_mle_1/sqrt(T);
% SE_22 = theta_mle_2 - 1.96*theta_mle_1/sqrt(T);
SEar1 = mean(SEAR1);
SEar_21 = theta_mle_AR1_2 + SEar1;
SEar_22 = theta_mle_AR1_2 - SEar1;
display(SEar_21);
display(SEar_22);

f1 = figure;
histfit(theta_mle_AR1(2,:),25,'kernel');
line([theta_mle_AR1_2, theta_mle_AR1_2], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([SEar_21 SEar_21 NaN SEar_22 SEar_22] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')

theta_mle_MA1_1 = mean(theta_mle_MA1(1,:));
% theta_mle_MA1_2 = mean((2*normcdf(theta_mle_MA1(2,:))-1));
theta_mle_MA1_2 = mean(theta_mle_MA1(2,:));
loglikeMA1 = mean(dLogLikMA1);
display(theta_mle_MA1_1);
display(theta_mle_MA1_2);

% SE_21 = theta_mle_2 + 1.96*theta_mle_1/sqrt(T);
% SE_22 = theta_mle_2 - 1.96*theta_mle_1/sqrt(T);
SEma1 = mean(SEMA1);
SEma_21 = theta_mle_MA1_2 + SEma1;
SEma_22 = theta_mle_MA1_2 - SEma1;
display(SEma_21);
display(SEma_22);

f2 = figure;
histfit(theta_mle_MA1(2,:),25,'kernel');
line([theta_mle_MA1_2, theta_mle_MA1_2], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([SEma_21 SEma_21 NaN SEma_22 SEma_22] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')

theta_mle_ARMA11_1 = mean(theta_mle_ARMA11(1,:));
% theta_mle_2 = mean((2*normcdf(theta_mle_AR1(2,:))-1));
theta_mle_ARMA11_2 = mean(theta_mle_ARMA11(2,:));
theta_mle_ARMA11_th = mean(theta_mle_ARMA11(3,:));
loglikeARMA11 = mean(dLogLikARMA11);
display(theta_mle_ARMA11_1);
display(theta_mle_ARMA11_2);
display(theta_mle_ARMA11_th);

% SE_21 = theta_mle_2 + 1.96*theta_mle_1/sqrt(T);
% SE_22 = theta_mle_2 - 1.96*theta_mle_1/sqrt(T);
SEarma11 = mean(SEARMA11);
SEarma11th = mean(SEARMA11th);
SEarma_21 = theta_mle_ARMA11_2 + SEarma11;
SEarma_22 = theta_mle_ARMA11_2 - SEarma11;
SEarma_21th = theta_mle_ARMA11_th + SEarma11th;
SEarma_22th = theta_mle_ARMA11_th - SEarma11th;
display(SEarma_21);
display(SEarma_22);
display(SEarma_21th);
display(SEarma_22th);

f3 = figure;
histfit(theta_mle_ARMA11(2,:),25,'kernel');
line([theta_mle_ARMA11_2, theta_mle_ARMA11_2], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([SEarma_21 SEarma_21 NaN SEarma_22 SEarma_22] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')
f4 = figure;
histfit(theta_mle_ARMA11(3,:),25,'kernel');
line([theta_mle_ARMA11_th, theta_mle_ARMA11_th], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([SEarma_21th SEarma_21th NaN SEarma_22th SEarma_22th] , [ylim NaN   ylim],'LineWidth', 0.5, 'Color', 'g','Displayname','St. Dev.')

%% Criteria Testing
 
%     Y_AIC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 1);
%     Y_BIC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 2);
%     Y_AICC(1) = loglikeAR1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 3);
%     
%     Y_AIC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 1);
%     Y_BIC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 2);
%     Y_AICC(2) = loglikeMA1(Y(:,1), [theta_mle_1 ; theta_mle_2], T, 3);    

    Y_AIC(1) = -2*loglikeAR1 + 4;
    Y_BIC(1) = -2*loglikeAR1 + 4*T/(T-3);
    Y_AICC(1) = -2*loglikeAR1 + 2*log(T)/T;
    
    Y_AIC(2) = -2*loglikeMA1 + 4;
    Y_BIC(2) = -2*loglikeMA1 + 4*T/(T-3);
    Y_AICC(2) = -2*loglikeMA1 + 2*log(T)/T; 
    
    Y_AIC(3) = -2*loglikeARMA11 + 6;
    Y_BIC(3) = -2*loglikeARMA11 + 6*T/(T-3);
    Y_AICC(3) = -2*loglikeARMA11 + 4*log(T)/T;

% end