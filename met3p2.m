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

ty = dates_(~(ret_ > ret_mean+4*ret_sd | ret_ < ret_mean-4*ret_sd));
n = size(Y,1);

% defining Fourier Freq 
omega = NaN(n,1);
meanY = mean(Y);

for k = (n-1)/2:-1:1
    omega((n-1)/2-k+1,1) = -2*pi()*k/n;
    omega((n-1)/2+k+1,1) = 2*pi()*k/n;    
end
omega((n-1)/2+1,1) = 0;

% finding all population autocov values
gammaY = zeros(n,1);
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
set(f12,'Visible','off');
autocorr(Y)
title('Returns: Autocorrelation')

f13 = figure;
set(f13,'Visible','off');
parcorr(Y)
title('Returns: Partial Autocorrelation')

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

%% Algorithm 

disp('testing algorithm')
estim(Y, gammaY);

%% 2. Models and Estimation

disp('Section 2');


%% 3. sdfsdf

disp('Section 3');


%% 4 Pockets of Predictability

disp('Section 4');


%%
%Save all plots
saveas(f11,'Figure 1.1.jpeg');
saveas(f12,'Figure 1.2.jpeg');
saveas(f13,'Figure 1.3.jpeg');
saveas(f14,'Figure 1.4.jpeg');
saveas(f15,'Figure 1.5.jpeg');

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

profile viewer
%profsave;
