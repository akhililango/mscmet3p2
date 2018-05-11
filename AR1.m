rng('default')     

% Monte-Carlo runs 
runs = 1000;

%AR(1) model specs
T = 300;
Y   = NaN(T,runs);
theta = [0 ; 0.9 ; 0.2]; %[mu_y phi_1 sigma]
p = 1;
q = 0;

%residuals
epsY = theta(3)*randn(T,runs);

%Generate the AR(1) process
Y(1,:) = epsY(1,:);
for t = 1:T-1
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
thetaStart = [0.1 ; 0.5 ; 0.15]; 
options = optimset('TolX', 0.0001, 'Display', 'off', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

%% DL/Inn
for i = 1:runs
%      objfun = @(thetaStart)(-loglikeARMA(thetaStart, Y(:,i), epsY(:,i), p, q));
     [YhatDL(:,i), meanY, GammaY] = estim(Y(:,i), thetaStart(p+q+2));
%      [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
%     disp(i)
end

%% MLE
r = sum((Y - YhatDL).^2,2);
r(1) = 1;
for i = 1:runs
%     objfun = @(thetaStart)(-loglikeARMA(thetaStart, Y(:,i), epsY, p, q));
    objfunTemp = 0;
    for m = 1:T
        objfunTemp = objfunTemp + ((Y(m,i) - YhatDL(m,i))^2)/r(m);
    end
    objfun = @(thetaStart)(log(((2*pi()*thetaStart(p+q+2))^(-T/2))*(prod(r)^(-1/2))*exp((-1/(2*(thetaStart(p+q+2))^2))*objfunTemp)));
        
    [theta_mle(:,i), dLogLik] = fminunc(objfun, thetaStart, options);
    disp(i)
end

%% Display

fprintf ('Log Likelihood value = %g \r', -dLogLik*runs);
display(mean(theta_mle(1,:)));
display(mean(theta_mle(2,:)));
display(mean(theta_mle(3,:)));
histfit(theta_mle(2,:));
