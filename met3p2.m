clear all
close all
clc
profile on

disp('Empirical pset 2, Metrics 3');

%% 0. Importing Data

data = readtable("AmazonData.xlsx");

ret_ = table2array(data(:,3));
dates_ = table2array(data(:,2));

ret_mean = mean(ret_);
ret_var = var(ret_);
ret_sd = sqrt(ret_var);

% Truncating the data (different from what the team did)
up_bound = ones(5273,1);
up_bound = ret_mean+4*ret_sd*up_bound;
low_bound = ones(5273,1);
low_bound = ret_mean-4*ret_sd*low_bound;
Y = min(ret_, up_bound);
Y = max(Y, low_bound);

ty = dates_;


% defining Fourier Freq 
n = size(Y,1);
omega = NaN(n,1);
meanY = mean(Y);

for k = (n-1)/2:-1:1
    omega((n-1)/2-k+1,1) = -2*pi()*k/n;
    omega((n-1)/2+k+1,1) = 2*pi()*k/n;    
end
omega((n-1)/2+1,1) = 0;

% finding all population autocov values
gammaY = zeros(n,1); %First element is actually what we call gammaY(0), and last one is gammaY(n-1)
for k = 1:n
    for t = 1:n+1-k
        gammaY(k) = gammaY(k) + (Y(t)-meanY)*(Y(t+k-1)-meanY);
    end
gammaY(k) = gammaY(k)/(n+1-k);
end
rhoY = gammaY/gammaY(1);

rhoYtest = autocorr(Y,n-1);
gammaYtest = rhoYtest*var(Y);

%% 1. Data Exploration

disp('Section 1');

f11 = figure;
set(f11,'Visible','off');
plot(ty,Y);
ax = gca;
ax.XTick = ty(1:200:end);
ax.XTickLabelRotation = 90;
datetick('x','yyyy','keepticks')
xlabel('Year');
ylabel('Returns');
title('Returns: Time series')
axis tight;
recessionplot;
saveas(f11,'Figure 1.1.jpeg');

f12 = figure;
set(f12,'Visible','off');
autocorr(Y)
title('Returns: Autocorrelation')
saveas(f12,'Figure 1.2.jpeg');

f13 = figure;
set(f13,'Visible','off');
parcorr(Y)
title('Returns: Partial Autocorrelation')
saveas(f13,'Figure 1.3.jpeg');

f14 = figure;
set(f14,'Visible','off');
histfit(Y,25,'kernel');
line([mean(Y), mean(Y)], ylim, 'LineWidth',1,'Color','r','LineStyle','-.')
line ([mean(Y)+std(Y) mean(Y)+std(Y) NaN mean(Y)-std(Y) ...
    mean(Y)-std(Y)] , [ylim NaN   ylim],'LineWidth', 0.5, ...
    'Color', 'g','Displayname','St. Dev.')
a=annotation('textbox',...
    [0.20 0.8 0.5 0.04],...
    'String',{'Mean = 0.0806', 'Variance = 11.3179', 'Skewness = 0.1421','Kurtosis = 6.6538'},...
    'FitBoxToText','on','LineStyle','none');
a.FontSize=10;
tf14=title('Returns: Empirical Density');
tf14.FontSize=10;
tf14.FontWeight = 'bold';
saveas(f14,'Figure 1.4.jpeg');



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
%% Periodogram from lecture
perio = zeros(n,1);

K = size(omega,1);
% for k = 1:K
%     for t = 1:n
%     eperio(:,n) = exp(1i*t*omega(k));
%     end
% end

for k= 1:K
    for t = 1:n
        perio(k) = perio(k) + Y(t)*(cos(omega(k)*t) - i*sin(omega(k)*t)); %exp(-1i*t*omega(k))
    end
    perio(k) = abs(perio(k))/n;
end

% plot
% xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})

f15 = figure
% set(f15,'Visible','off');
plot(omega,perio)
ax = gca;
ax.XTick = [-pi -pi/2 -pi/4 0 pi/4 pi/2 pi]
ax.XTickLabel = {'-\pi -\pi/2 -\pi/4 0 \pi/4 \pi/2 \pi'}
ax.XTickLabelMode = 'auto'
% set(gca,'XTicks',omega); % Change x-axis ticks
% set(gca,'XTickLabels',omega); % Change x-axis ticks labels to desired values.
xlabel('Frequency')
ylabel('Density')
title('Periodogram')
saveas(f15,'Figure 1.5.jpeg');

ftest = figure
% set(f15,'Visible','off');
plot(omega,perio)


disp('graphs created')


%% 2. Models and Estimation

disp('Section 2');

Y_sim = NaN(1000,1);
Y_sim(1) = 0;
epsilon = normrnd(0,0.2,[1000,1]);
epsilon(1) = NaN;
for t = 1:999
    Y_sim(t+1)=0.9*Y_sim(t)+epsilon(t+1);
end

f_ar1=figure
plot(Y_sim)
set(f_ar1, 'Visible','on')
