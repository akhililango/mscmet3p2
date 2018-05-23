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

f15 = figure
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

    YhatWN(t+1) = theta_WN_2 + theta_WN_1*randn();
    
    MSFEWN(t-99) = sum((Y(101:t+1) - YhatWN(101:t+1)).^2)/(t-99);
    MAFEWN(t-99) = sum(abs(Y(101:t+1) - YhatWN(101:t+1)))/(t-99);
end

%% AR(1)

YhatAR1 = Y(1:100);
for t = 100:T-1
    thetaStart = [0.1 ; 0.5]; 
    options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
    objfun = @(thetaStart)(-loglikeAR1(Y(t-99:t), thetaStart, 100));
    [theta_AR1, dLogLik] = fminunc(objfun, thetaStart, options);

    theta_AR1_1 = theta_AR1(1);
    theta_AR1_2 = theta_AR1(2);
    % display(theta_AR1_1);
    % display(theta_AR1_2);
    % 
    % SE_AR1_21 = theta_AR1_2 + 1.96*theta_AR1_1/sqrt(100);
% SE_AR1_22 = theta_AR1_2 - 1.96*theta_AR1_1/sqrt(100);

    Y_AIC(1,t-99) = loglikeAR1(Y(t-99:t), [theta_AR1_1 ; theta_AR1_2], 100, 1);
    Y_BIC(1,t-99) = loglikeAR1(Y(t-99:t), [theta_AR1_1 ; theta_AR1_2], 100, 2);
    Y_AICC(1,t-99) = loglikeAR1(Y(t-99:t), [theta_AR1_1 ; theta_AR1_2], 100, 3);

    YhatAR1(t+1) = theta_AR1_2*YhatAR1(t);
    
     MSFEAR1(t-99) = sum((Y(101:t+1) - YhatAR1(101:t+1)).^2)/(t-99);  % or /T
     MAFEAR1(t-99) = sum(abs(Y(101:t+1) - YhatAR1(101:t+1)))/(t-99);  % or /T
end
 
% dmAR1.stat = dm.test(MSFEAR1,MAFEAR1)$statistic
% dmAR1.pval = dm.test(MSFEAR1,MAFEAR1)$p.value

%% MA(1)

    YhatMA1 = Y(1:100);
for t = 100:T-1
    thetaStart = [0.1 ; 0.5]; 
    options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
    objfun = @(thetaStart)(-loglikeMA1(Y(t-99:t), thetaStart, 100));
    [theta_MA1, dLogLik] = fminunc(objfun, thetaStart, options);

    theta_MA1_1 = theta_MA1(1);
    theta_MA1_2 = theta_MA1(2);
    % display(theta_MA1_1);
    % display(theta_MA1_2);
    % 
    % SE_MA1_21 = theta_MA1_2 + 1.96*theta_MA1_1/sqrt(100);
    % SE_MA1_22 = theta_MA1_2 - 1.96*theta_MA1_1/sqrt(100);

    Y_AIC(2,t-99) = loglikeMA1(Y(t-99:t), [theta_MA1_1 ; theta_MA1_2], 100, 1);
    Y_BIC(2,t-99) = loglikeMA1(Y(t-99:t), [theta_MA1_1 ; theta_MA1_2], 100, 2);
    Y_AICC(2,t-99) = loglikeMA1(Y(t-99:t), [theta_MA1_1 ; theta_MA1_2], 100, 3);
    

    YhatMA1(t+1) = theta_MA1_2*(theta_MA1_2 + theta_MA1_1*randn());
    
    MSFEMA1(t-99) = sum((Y(101:t+1) - YhatMA1(101:t+1)).^2)/(t-99);  % or /T
    MAFEMA1(t-99) = sum(abs(Y(101:t+1) - YhatMA1(101:t+1)))/(t-99);  % or /T
end

    
%% ARMA(1,1)

YhatARMA11 = Y(1:100);
for t = 100:T-1
    thetaStart = [0.1 ; 0.4 ; -0.5]; 
    options = optimset('TolX', 0.0001, 'Display', 'iter-detailed', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
    objfun = @(thetaStart)(-loglikeARMA11(Y(t-99:t), thetaStart, 100));
    [theta_ARMA11, dLogLik] = fminunc(objfun, thetaStart, options);

    theta_ARMA11_1 = theta_ARMA11(1);
    theta_ARMA11_2 = theta_ARMA11(2);
    theta_ARMA11_3 = theta_ARMA11(3);
    % display(theta_ARMA11_1);
    % display(theta_ARMA11_2);
    % display(theta_ARMA11_3);
    % 
    % SE_ARMA11_21 = theta_ARMA11_2 + 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_22 = theta_ARMA11_2 - 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_31 = theta_ARMA11_3 + 1.96*theta_ARMA11_1/sqrt(100);
    % SE_ARMA11_32 = theta_ARMA11_3 - 1.96*theta_ARMA11_1/sqrt(100);


    Y_AIC(3,t-99) = loglikeARMA11(Y(t-99:t), [theta_ARMA11_1 ; theta_ARMA11_2 ; theta_ARMA11_3], 100, 1);
    Y_BIC(3,t-99) = loglikeARMA11(Y(t-99:t), [theta_ARMA11_1 ; theta_ARMA11_2 ; theta_ARMA11_3], 100, 2);
    Y_AICC(3,t-99) = loglikeARMA11(Y(t-99:t), [theta_ARMA11_1 ; theta_ARMA11_2 ; theta_ARMA11_3], 100, 3);
    
    YhatARMA11(t+1) = theta_ARMA11_2*YhatARMA11(t) + theta_ARMA11_3*(theta_ARMA11_3 + theta_ARMA11_1*randn()); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ???????????
    
    MSFEARMA11(t-99) = sum((Y(101:t+1) - YhatARMA11(101:t+1)).^2)/(t-99);  % or /T
    MAFEARMA11(t-99) = sum(abs(Y(101:t+1) - YhatARMA11(101:t+1)))/(t-99);  % or /T
end

  
%% DM Test

[dm_WN_stat, dm_WN_pval] = dieboldmariano(MSFEWN,MAFEWN);
[dm_AR1_stat, dm_AR1_pval] = dieboldmariano(MSFEAR1,MAFEAR1);
[dm_MA1_stat, dm_MA1_pval] = dieboldmariano(MSFEMA1,MAFEMA1);
[dm_ARMA11_stat, dm_ARMA11_pval] = dieboldmariano(MSFEARMA11,MAFEARMA11);
  
%% 4 Pockets of Predictability

disp('Section 4');

rl_AR1 = YhatAR1(101:T) - YhatWN(101:T);
rl_MA1 = YhatMA1(101:T) - YhatWN(101:T);
% rl_ARMA11 = YhatARMA11(101:T) - YhatWN(101:T);

  lsAR1 = zeros(1,T-301);
  lsMA1 = zeros(1,T-301);
%   lsARMA11 = zeros(1,T-300);

for s = 1:T-301
  lsAR1(1,s) = sum(rl_AR1(s:s+201));
  lsMA1(1,s) = sum(rl_MA1(s:s+201));
%   lsARMA11(s-200) = sum(lARMA11((s-100):(s+100)));
end

  xAR1 = zeros(1,T-301);
  xAR1(lsAR1>0) = 1;
  xAR1(lsAR1<0) = -1;
  xMA1 = zeros(1,T-301);
  xMA1(lsMA1>0) = 1;
  xMA1(lsMA1<0) = -1;
%   xARMA11 = ones(1,T-300);
%   xARMA11(lsARMA11<0) = -1;

f41 = figure;
plot(xAR1,lsAR1);
f42 = figure;
plot(xMA1,lsMA1);
% f43 = figure;
% plot(xARMA11,lsARMA11);

%% DM Test

[dms_AR1_stat, dms_AR1_pval] = dieboldmariano(MSFEAR1((lsAR1 < -200|lsAR1 > 200)),MAFEAR1((lsAR1 < -200|lsAR1 > 200)));
[dms_MA1_stat, dms_MA1_pval] = dieboldmariano(MSFEMA1((lsMA1 < -200|lsMA1 > 200)),MAFEMA1((lsMA1 < -200|lsMA1 > 200)));
% [dm_ARMA11_stat, dm_ARMA11_pval] = dieboldmariano(MSFEARMA11,MAFEARMA11);

disp('No, they didnt ouperform. The conditional mean model is hard to beat in out-of-sample forecasating.');
disp('The DM test should never be used in this way because it is for Asymptotic testing and this has not been corrected for finite samples');
%%
%Save all plots
saveas(f11,'Figure 1.1.jpeg');
saveas(f12,'Figure 1.2.jpeg');
saveas(f13,'Figure 1.3.jpeg');
saveas(f14,'Figure 1.4.jpeg');
saveas(f15,'Figure 1.5.jpeg');
saveas(f41,'Figure 4.1.jpeg');
saveas(f42,'Figure 4.2.jpeg');
% saveas(f43,'Figure 4.3.jpeg');

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