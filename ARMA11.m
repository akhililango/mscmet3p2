clear all
close all
clc
disp('ARMA1')
rng('default')  

% Monte-Carlo runs 
runs = 1000;

%ARMA(1,1) model specs
T = 300;
Y   = NaN(T,runs);
theta = [0 ; 0.4 ; 0.3; -0.2]; 
p = 1;
q = 1;

%residuals
epsY = theta(2)*randn(T,runs);

%Generate the ARMA(1,1) process
Y(1,:) = epsY(1);
for t = 1:T-1
   Y(t+1,:) = theta(1) + theta(3)*Y(t,:) + theta(4)*epsY(t,:) + epsY(t+1,:);
   
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
thetaStart = [0 ; 0.5 ; 0.2 ; -0.15]; rng('default')   
options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

for i = 1:runs
    objfun = @(thetaStart)(-loglikeARMA11(Y(:,i)-mean(Y(:,i)), thetaStart, T));
    [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
end

fprintf ('Log Likelihood value = %g \r', -dLogLik*runs);
% display(mean(theta_mle(1,:)));
display(mean(theta_mle(2,:)));
display(mean(theta_mle(3,:)));
display(mean(theta_mle(4,:)));
histfit(theta_mle(2,:));
histfit(theta_mle(3,:));
  