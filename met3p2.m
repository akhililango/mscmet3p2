clear all
close all
clc
% profile on

disp('Empirical pset 2, Metrics 3');

%% 0. Importing Data

data = readtable("AmazonData.xlsx");

ret_ = table2array(data(:,3));
dates_ = table2array(data(:,2));

ret_mean = mean(ret_);
ret_var = var(ret_);
ret_sd = sqrt(ret_var);
Y = ret_;
Y(ret_ > ret_mean+4*ret_sd) = ret_mean+4*ret_sd;
Y(ret_ < ret_mean-4*ret_sd) = ret_mean-4*ret_sd;
ty = dates_;
T = size(Y,1);

% defining Fourier Freq 
omega = zeros(T,1);
for k = (T-1)/2:-1:1
    omega((T-1)/2-k+1,1) = -2*pi()*k/T;
    omega((T-1)/2+k+1,1) = 2*pi()*k/T;    
end
omega((T-1)/2+1,1) = 0;

% finding all sample autocov values
meanY = mean(Y);
gammaY = zeros(T,1);
for k = 1:T
    for t = 1:T+1-k
        gammaY(k) = gammaY(k) + (Y(t)-meanY)*(Y(t+k-1)-meanY);
    end
    gammaY(k) = gammaY(k)/(T+1-k);
end
rhoY = gammaY/gammaY(1); %sample autocorrelation

rhoYtest = autocorr(Y,T-1);
gammaYtest = rhoYtest*var(Y);

%% 1. Data Exploration

disp('Section 1');

f11 = figure;
%set(f11,'Visible','off');
plot(ty,Y,'LineWidth',2);
ax = gca;
ax.XTick = ty(1:200:end);
ax.XTickLabelRotation = 90;
datetick('x','yyyy','keepticks')
xlabel('Year');
ylabel('Returns');
title('Returns: Time series')
axis tight;
recessionplot;

f12 = figure;
%set(f12,'Visible','off');
autocorr(Y)
title('Returns: Autocorrelation')

f13 = figure;
%set(f13,'Visible','off');
parcorr(Y)
title('Returns: Partial Autocorrelation')

f14 = figure;
%set(f14,'Visible','off');
histfit(Y,25,'kernel')
line([mean(Y), mean(Y)], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([mean(Y)+std(Y) mean(Y)+std(Y) NaN mean(Y)-std(Y) ...
    mean(Y)-std(Y)] , [ylim NaN   ylim],'LineWidth', 0.5, ...
    'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.20 0.8 0.5 0.04],...
    'String',{'Mean = 0.0806', 'Variance = 11.3179', 'Skewness = 0.1421','Kurtosis = 6.6538'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=10;
ty=title('Returns: Empirical Density');
ty.FontSize=10;
ty.FontWeight = 'bold';


% % % %  %
% f15 = figure;
% set(f15,'Visible','off');
% periodogram(ret);
% title('Returns: Periodogram')

% datevar = datenum(t);
% f15 = figure;
% subplot(2,1,1); 
% plot(datevar,Y); 
% subplot(2,1,2); 
% PowerSpectrum=PlotFrequencySpectrum(datevar,Y,3,1,0);
% hold on
%% 1.8. Periodogram from lecture
perio = zeros(T,1);

K = size(omega,1);
% for k = 1:K
%     for t = 1:n
%     eperio(:,n) = exp(1i*t*omega(k));
%     end
% end

for k= 1:K
    for t = 1:T
        perio(k) = perio(k) + Y(t)*(cos(omega(k)*t) - i*sin(omega(k)*t)); %exp(-1i*t*omega(k))
    end
    perio(k) = abs(perio(k))/T;
end

% plot
% xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})

f15 = figure;
% set(f15,'Visible','off');
plot(omega,perio,'LineWidth',2)
ax = gca;
ax.XTick = [-pi -pi/2 -pi/4 0 pi/4 pi/2 pi]
ax.XTickLabel = {'-\pi -\pi/2 -\pi/4 0 \pi/4 \pi/2 \pi'}
ax.XTickLabelMode = 'auto'
% set(gca,'XTicks',omega); % Change x-axis ticks
% set(gca,'XTickLabels',omega); % Change x-axis ticks labels to desired values.
xlabel('Frequency')
ylabel('Density')
title('Periodogram')

disp('graphs created')


%% 2. Models and Estimation

disp('Section 2');


%% 3. Forecasting

disp('Section 3');

%% WN

YhatWN = Y(1:100);
for t = 100:T-1
    theta_WN_1 = std(Y(t-99:t));
    theta_WN_2 = mean(Y(t-99:t));
%     display(theta_WN1_1);
%     display(theta_WN1_2);

%     SE_WN_1 = theta_WN1_2 + 1.96*theta_WN1_1/sqrt(100);
%     SE_WN_2 = theta_WN1_2 - 1.96*theta_WN1_1/sqrt(100);

    YhatWN(t+1) = theta_WN_2;
end

%     MSFEWN(t-99) = sum((Y(101:t+1) - YhatWN(101:t+1)).^2)/(t-99);
%     MAFEWN(t-99) = sum(abs(Y(101:t+1) - YhatWN(101:t+1)))/(t-99);

    MSFEWN = sum((Y(101:T) - YhatWN(101:T)).^2)/(T-100);
    MAFEWN = sum(abs(Y(101:T) - YhatWN(101:T)))/(T-100);
    dWNF = (Y(101:T) - YhatWN(101:T)).^2;
    dWNA = abs(Y(101:T) - YhatWN(101:T));
    
%% AR(1)

YhatAR1 = Y(1:100);
for t = 100:T-1
%     thetaStart = [0.1 ; 0.5]; 
%     options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
%     objfun = @(thetaStart)(-loglikeAR1(Y(t-99:t)-mean(Y(t-99:t)), thetaStart, 100));
%     [theta_AR1, dLogLik] = fminunc(objfun, thetaStart, options);
% 
%     theta_AR1_1 = theta_AR1(1);
%     theta_AR1_2 = theta_AR1(2);
%     loglike = dLogLik;
    % display(theta_AR1_1);
    % display(theta_AR1_2);
    % 
    % SE_AR1_21 = theta_AR1_2 + 1.96*theta_AR1_1/sqrt(100);
%     SE_AR1_22 = theta_AR1_2 - 1.96*theta_AR1_1/sqrt(100);

        
    ToEstMdl = arima(1,0,0);
    [EstMdl,EstParamCov,logL,info] = estimate(ToEstMdl,Y(t-99:t)-mean(Y(t-99:t)));
    theta_AR1_1 = sqrt(info.X(3));
    theta_AR1_2 = info.X(2);
    loglike = logL;

    Y_AIC(1,t-99) = -2*loglike + 4;
    Y_BIC(1,t-99) = -2*loglike + 4*T/(T-3);
    Y_AICC(1,t-99) = -2*loglike + 2*log(T)/T;

    YhatAR1(t+1) = theta_AR1_2*YhatAR1(t);
end
 
%      MSFEAR1(t-99) = sum((Y(101:t+1) - YhatAR1(101:t+1)).^2)/(t-99);  % or /T
%      MAFEAR1(t-99) = sum(abs(Y(101:t+1) - YhatAR1(101:t+1)))/(t-99);  % or /T
     
     MSFEAR1 = sum((Y(101:T) - YhatAR1(101:T)).^2)/(T-100);  % or /T
     MAFEAR1 = sum(abs(Y(101:T) - YhatAR1(101:T)))/(T-100);  % or /T
     dAR1F = (Y(101:T) - YhatAR1(101:T)).^2 - (Y(101:T) - YhatWN(101:T)).^2;
     dAR1A = abs(Y(101:T) - YhatAR1(101:T)) - abs(Y(101:T) - YhatWN(101:T));
     Y_AIC_AR1 = mean(Y_AIC);
     Y_BIC_AR1 = mean(Y_BIC);
     Y_AICC_AR1 = mean(Y_AICC);

     
% dmAR1.stat = dm.test(MSFEAR1,MAFEAR1)$statistic
% dmAR1.pval = dm.test(MSFEAR1,MAFEAR1)$p.value

%% MA(1)

    YhatMA1 = Y(1:100);
for t = 100:T-1
%     thetaStart = [0.1 ; 0.5]; 
%     options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
%     objfun = @(thetaStart)(-loglikeMA1(Y(t-99:t)-mean(Y(t-99:t)), thetaStart, 100));
%     [theta_MA1, dLogLik] = fminunc(objfun, thetaStart, options);
% 
%     theta_MA1_1 = theta_MA1(1);
%     theta_MA1_2 = theta_MA1(2);
%     loglike = dLogLik;
    % display(theta_MA1_1);
    % display(theta_MA1_2);
    % 
    % SE_MA1_21 = theta_MA1_2 + 1.96*theta_MA1_1/sqrt(100);
    % SE_MA1_22 = theta_MA1_2 - 1.96*theta_MA1_1/sqrt(100);
    
    ToEstMdl = arima(0,0,1);
    [EstMdl,EstParamCov,logL,info] = estimate(ToEstMdl,Y(t-99:t)-mean(Y(t-99:t)));
    theta_MA1_1 = sqrt(info.X(3));
    theta_MA1_2 = info.X(2);
    loglike = logL;
    
    Y_AIC(1,t-99) = -2*loglike + 4;
    Y_BIC(1,t-99) = -2*loglike + 4*T/(T-3);
    Y_AICC(1,t-99) = -2*loglike + 2*log(T)/T;

    YhatMA1(t+1) = theta_MA1_2*(theta_MA1_2 + theta_MA1_1*randn());
    
%     MSFEMA1(t-99) = sum((Y(101:t+1) - YhatMA1(101:t+1)).^2)/(t-99);  % or /T
%     MAFEMA1(t-99) = sum(abs(Y(101:t+1) - YhatMA1(101:t+1)))/(t-99);  % or /T
end
    
        
    MSFEMA1 = sum((Y(101:T) - YhatMA1(101:T)).^2)/(T-100);  % or /T
    MAFEMA1 = sum(abs(Y(101:T) - YhatMA1(101:T)))/(T-100);  % or /T
    dMA1F = (Y(101:T) - YhatMA1(101:T)).^2 - (Y(101:T) - YhatWN(101:T)).^2;
    dMA1A = abs(Y(101:T) - YhatMA1(101:T)) - abs(Y(101:T) - YhatWN(101:T));
    Y_AIC_MA1 = mean(Y_AIC);
    Y_BIC_MA1 = mean(Y_BIC);
    Y_AICC_MA1 = mean(Y_AICC);
    
%% ARMA(1,1)

YhatARMA11 = Y(1:100);
for t = 100:T-1 
%     thetaStart = [0.1 ; 0.4 ; -0.5]; 
%     options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
%     objfun = @(thetaStart)(-loglikeARMA11(Y(t-99:t)-mean(Y(t-99:t)), thetaStart, 100));
%     [theta_ARMA11, dLogLikARMA11] = fminunc(objfun, thetaStart, options);
% 
%         
% %     res = Hevia_arma_mle(Y, 1, 1);
% %     theta_mle_ARMA11(1) = res.sigma;
% %     theta_mle_ARMA11(2) = 2*normcdf(res.ar)-1;
% %     theta_mle_ARMA11(3) = 2*normcdf(res.ma)-1;
% %     invhess = inv(res.hess);
% %     SEARMA11 = 1.96*invhess(1,1);
% %     SEARMA11th = 1.96*invhess(2,2);
% %     dLogLikARMA11 = res.loglike;
% %     
%     theta_ARMA11_1 = theta_ARMA11(1);
%     theta_ARMA11_2 = theta_ARMA11(2);
%     theta_ARMA11_3 = theta_ARMA11(3);
%     loglike = dLogLikARMA11;
    % display(theta_ARMA11_1);
    % display(theta_ARMA11_2);
    % display(theta_ARMA11_3);
    % 
    % SE_ARMA11_21 = theta_ARMA11_2 + 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_22 = theta_ARMA11_2 - 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_31 = theta_ARMA11_3 + 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_32 = theta_ARMA11_3 - 1.96*theta_ARMA11_1/sqrt(100);
    
    
    ToEstMdl = arima(1,0,1);
    [EstMdl,EstParamCov,logL,info] = estimate(ToEstMdl,Y(t-99:t)-mean(Y(t-99:t)));
    if ~(size(info.X,1)<4)
    theta_ARMA11_1 = sqrt(info.X(4));
    theta_ARMA11_2 = info.X(2);
    theta_ARMA11_3 = info.X(3);
    loglike = logL;
    end
    Y_AIC(3) = 2*loglike + 6;
    Y_BIC(3) = 2*loglike + 6*T/(T-3);
    Y_AICC(3) = 2*loglike + 4*log(T)/T;

    YhatARMA11(t+1) = theta_ARMA11_2*YhatARMA11(t) + theta_ARMA11_3*(theta_ARMA11_3 + theta_ARMA11_1*randn()); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ???????????

%     MSFEARMA11(t-99) = sum((Y(101:t+1) - YhatARMA11(101:t+1)).^2)/(t-99);  % or /T
%     MAFEARMA11(t-99) = sum(abs(Y(101:t+1) - YhatARMA11(101:t+1)))/(t-99);  % or /T
end
    
    MSFEARMA11 = sum((Y(101:T) - YhatARMA11(101:T)).^2)/(T-100);  % or /T
    MAFEARMA11 = sum(abs(Y(101:T) - YhatARMA11(101:T)))/(T-100);  % or /T
    dARMA11F = (Y(101:T) - YhatARMA11(101:T)).^2 - (Y(101:T) - YhatWN(101:T)).^2;
    dARMA11A = abs(Y(101:T) - YhatARMA11(101:T)) - abs(Y(101:T) - YhatWN(101:T));
    Y_AIC_ARMA11 = mean(Y_AIC);
    Y_BIC_ARMA11 = mean(Y_BIC);
    Y_AICC_ARMA11 = mean(Y_AICC);
  
%% DM Test

% defining Fourier Freq 
dT = T - 100;
omegadm = zeros(dT,1);
for k = (dT-1)/2:-1:1
    omegadm((dT-1)/2-k+1,1) = -2*pi()*k/dT;
    omegadm((dT-1)/2+k+1,1) = 2*pi()*k/dT;    
end
omegadm((dT-1)/2+1,1) = 0;

% periodm = zeros(dT,1);
% 
% K = size(omegadm,1);
% % for k = 1:K
% %     for t = 1:n
% %     eperio(:,n) = exp(1i*t*omega(k));
% %     end
% % end
% 
% for k= 1:K
%     for t = 1:dT
%         periodm(k) = periodm(k) + dAR1F(t)*(cos(omegadm(k)*t) - i*sin(omegadm(k)*t)); %exp(-1i*t*omega(k))
%     end
%     periodm(k) = abs(periodm(k))/dT;
% end

mid = (dT-1)/2+1;

[dm_AR1F_stat, dm_AR1F_pval] = dieboldmariano(dAR1F,omegadm,mid);
[dm_MA1F_stat, dm_MA1F_pval] = dieboldmariano(dMA1F,omegadm,mid);
[dm_ARMA11F_stat, dm_ARMA11F_pval] = dieboldmariano(dARMA11F,omegadm,mid);

[dm_AR1A_stat, dm_AR1A_pval] = dieboldmariano(dAR1A,omegadm,mid);
[dm_MA1A_stat, dm_MA1A_pval] = dieboldmariano(dMA1A,omegadm,mid);
[dm_ARMA11A_stat, dm_ARMA11A_pval] = dieboldmariano(dARMA11A,omegadm,mid);

%% 4 Pockets of Predictability

disp('Section 4');

rl_AR1 = (YhatAR1(101:T) - Y(101:T)).^2 - (YhatWN(101:T) - Y(101:T)).^2;
rl_MA1 = (YhatMA1(101:T) - Y(101:T)).^2 - (YhatWN(101:T) - Y(101:T)).^2;
rl_ARMA11 = (YhatARMA11(101:T) - Y(101:T)).^2 - (YhatWN(101:T) - Y(101:T)).^2;

% rl_AR1 = (YhatAR1(101:T) - YhatWN(101:T)).^2 - (YhatWN(101:T) - YhatWN(101:T)).^2;
% rl_MA1 = (YhatMA1(101:T) - YhatWN(101:T)).^2 - (YhatWN(101:T) - YhatWN(101:T)).^2;
% rl_ARMA11 = (YhatARMA11(101:T) - YhatWN(101:T)).^2 - (YhatWN(101:T) - YhatWN(101:T)).^2;

  lsAR1 = zeros(1,T-301);
  lsMA1 = zeros(1,T-301);
  lsARMA11 = zeros(1,T-301);

for s = 1:T-301
  lsAR1(1,s) = sum(rl_AR1(s:s+201));
  lsMA1(1,s) = sum(rl_MA1(s:s+201));
  lsARMA11(1,s) = sum(rl_ARMA11(s:s+201));
end

  xAR1 = zeros(1,T-301);
  xAR1(lsAR1>0) = 1;
  xAR1(lsAR1<0) = -1;
  
  xMA1 = zeros(1,T-301);
  xMA1(lsMA1>0) = 1;
  xMA1(lsMA1<0) = -1;
  
  xARMA11 = ones(1,T-301);
  xARMA11(lsARMA11>0) = 1;
  xARMA11(lsARMA11<0) = -1;

f41 = figure;
plot(1:T-301,xAR1,1:T-301,lsAR1);
f42 = figure;
plot(1:T-301,xAR1,1:T-301,lsMA1);
f43 = figure;
plot(1:T-301,xARMA11,1:T-301,lsARMA11);

%%

%%%%% AR1
ind_AR1(1) = 1;
i = 2;
for t = 1:T-302
    if lsAR1(t)*lsAR1(t+1)<0
        ind_AR1(i) = t;
        i = i+1;
    end
end
ind_AR1(i) = T-301;

t0 = 1;
sumlsAR1 = zeros(size(ind_AR1,2)-1,1);
for i = 1:size(ind_AR1,2)-1
    sumlsAR1(i) = 0;
    for t = t0:ind_AR1(i+1)
        sumlsAR1(i) = sumlsAR1(i) + lsAR1(t);
    end
    t0 = ind_AR1(i+1)+1;
end
% sumlsAR1 = abs(sumlsAR1);
[~,ilsAR1] = min(sumlsAR1);

%%%%%% MA1
ind_MA1(1) = 1;
i = 2;
for t = 1:T-302
    if lsMA1(t)*lsMA1(t+1)<0
        ind_MA1(i) = t;
        i = i+1;
    end
end
ind_MA1(i) = T-301;

t0 = 1;
sumlsMA1 = zeros(size(ind_MA1,2)-1,1);
for i = 1:size(ind_MA1,2)-1
    sumlsMA1(i) = 0;
    for t = t0:ind_MA1(i+1)
        sumlsMA1(i) = sumlsMA1(i) + lsAR1(t);
    end
    t0 = ind_MA1(i+1)+1;
end
% sumlsMA1 = abs(sumlsMA1);
[~,ilsMA1] = min(sumlsMA1);

%%%%%% ARMA11
ind_ARMA11(1) = 1;
i = 2;
for t = 1:T-302
    if lsARMA11(t)*lsARMA11(t+1)<0
        ind_ARMA11(i) = t;
        i = i+1;
    end
end
ind_ARMA11(i) = T-301;

t0 = 1;
sumlsARMA11 = zeros(size(ind_ARMA11,2)-1,1);
for i = 1:size(ind_ARMA11,2)-1
    sumlsARMA11(i) = 0;
    for t = t0:ind_ARMA11(i+1)
        sumlsARMA11(i) = sumlsARMA11(i) + lsAR1(t);
    end
    t0 = ind_ARMA11(i+1)+1;
end
sumlsARMA11 = abs(sumlsARMA11);
[~,ilsARMA11] = max(sumlsARMA11);

%% DM Test

%AR1
dT = ind_AR1(ilsAR1+1) - ind_AR1(ilsAR1);
omegadms_AR1 = zeros(dT,1);
for k = (dT-1)/2:-1:1
    omegadms_AR1((dT-1)/2-k+1,1) = -2*pi()*k/dT;
    omegadms_AR1((dT-1)/2+k+1,1) = 2*pi()*k/dT;    
end
omegadms_AR1((dT-1)/2+1,1) = 0;
mid = ceil((dT-1)/2+1);
[dms_AR1_stat, dms_AR1_pval] = dieboldmariano(lsAR1(ind_AR1(ilsAR1)+1:ind_AR1(ilsAR1+1)), omegadms_AR1, mid);

%MA1
dT = ind_MA1(ilsMA1+1) - ind_MA1(ilsMA1);
omegadms_MA1 = zeros(dT,1);
for k = (dT-1)/2:-1:1
    omegadms_MA1((dT-1)/2-k+1,1) = -2*pi()*k/dT;
    omegadms_MA1((dT-1)/2+k+1,1) = 2*pi()*k/dT;    
end
omegadms_MA1(ceil((dT-1)/2+1),1) = 0;
mid = ceil((dT-1)/2+1);
[dms_MA1_stat, dms_MA1_pval] = dieboldmariano(lsMA1(ind_MA1(ilsMA1)+1:ind_MA1(ilsMA1+1)), omegadms_MA1, mid);

%ARMA11
dT = ind_ARMA11(ilsARMA11+1) - ind_ARMA11(ilsARMA11);
for k = (dT-1)/2:-1:1
    omegadms_ARMA11((dT-1)/2-k+1,1) = -2*pi()*k/dT;
    omegadms_ARMA11((dT-1)/2+k+1,1) = 2*pi()*k/dT;    
end
omegadms_ARMA11(ceil((dT-1)/2+1),1) = 0;
mid = ceil((dT-1)/2+1);
[dms_ARMA11_stat, dms_ARMA11_pval] = dieboldmariano(lsARMA11(ind_ARMA11(ilsARMA11)+1:ind_ARMA11(ilsARMA11+1)), omegadms_ARMA11, mid);

disp('No, they didnt ouperform. The conditional mean model is hard to beat in out-of-sample forecasating.');
disp('The DM test should never be used in this way because it is for Asymptotic testing and has not been corrected for finite samples');
%%
%Save all plots
saveas(f11,'Figure 1.1.jpeg');
saveas(f12,'Figure 1.2.jpeg');
saveas(f13,'Figure 1.3.jpeg');
saveas(f14,'Figure 1.4.jpeg');
saveas(f15,'Figure 1.5.jpeg');
saveas(f41,'Figure 4.1.jpeg');
saveas(f42,'Figure 4.2.jpeg');
saveas(f43,'Figure 4.3.jpeg');

% saveas(f211,'Figure 2.1.1.png');
% saveas(f212,'Figure 2.1.2.png');
% saveas(f22,'Figure 2.2.png');
% saveas(f31,'Figure 3.1.png');
% saveas(f32,'Figure 3.2.png');
% saveas(f33,'Figure 3.3.png');
% saveas(f41,'Figure 4.1.png');
% saveas(f42,'Figure 4.2.png');
% saveas(f43,'Figure 4.3.png');
% saveas(f44,'Figure 4.4.png');

% profile viewer
%profsave;