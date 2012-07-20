% Filename: AfghanSimpleStats.m
% Author: Andrew Zammit Mangion
% Date: February 2011
% 
% Description: Performs some elementary stats on the temporal evolution of
% events in the AWD, including a Shapiro-Wilks test for normality, a
% Levene's test for homoscedasticity of the increments.
%
% Requires: STATISTICS TOOLBOX
%           ../AfghanDataAllDay (AWD in Matlab format)
%           swtest.m (Available from Matlab Exchange)
%           Levenetest.m (Available from Matlab Exchange)

clear all
% Load dataset
load('../AfghanDataAllDay')

%Sort into weeks
numofweeks = ceil(length(spikeAll)/7);
Weeknum = zeros(numofweeks,1);
shift = [57.5,28.5]; % Shift and scale to warp Afghanistan into roughly a 36 x 36 square
scale = [2,3.3];
i=1;
k=1;
% hold on
while i < length(spikeAll)
   spikeAllWeekData(k).Coords = [];
   for j = 0:6  %Cluster in weeks and put on 30x30 square
    if length(spikeAll(i).Coords) > 0
        spikeAll(i).Coords(:,1) = (spikeAll(i).Coords(:,1) - shift(1))*scale(1);
        spikeAll(i).Coords(:,2) = (spikeAll(i).Coords(:,2) - shift(2))*scale(2);
        spikeAllWeekData(k).Coords = [spikeAllWeekData(k).Coords;spikeAll(i).Coords];
    end
    Daynum(i) = size(spikeAll(i).Coords,1);
    Weeknum(k) = Weeknum(k) + length(spikeAll(i).Coords);
    i = i+1;
   end
   hold on
   k = k+1;
end

Weeknum(end) = []; % Incomplete final week

% Plot time series
figure('Position',[100,100,750,350])
plot(Weeknum,'k')
hold on 
weekinter = 365/7;
for i = 1:6
    yearstr = num2str(2003+i);
    text(round((i-1)*weekinter)+20,1100,yearstr);
    plot([round(i*weekinter),round(i*weekinter)],[0,1200],'k--')
end
xlabel('Week number','Fontsize',14)
ylabel('Number of activity reports','Fontsize',14)
axis([0, length(Weeknum), 0 1200])

set(gcf,'PaperPositionMode','auto')
print -dpng -r600 ../../Figures/AllGrowth.png

%Shapiro-Wilk test
for i = 1:length(Weeknum)-1
    Stat(i) = (Weeknum(i+1) - Weeknum(i))/Weeknum(i);
end
subplot(2,1,1)
histfit(Stat)
%Here we notice there are a lot of outliers - and we remove those noisy
%values. A strict outlier test might be adequate here.
i=1;
j=1;
indexremoved = [];
while i <= length(Stat)
    if(Stat(i) < -0.65) || (Stat(i) > 0.65)
        Stat(i) = [];
        indexremoved(end+1) = j;
    else
        i = i+1;
    end
    j= j+1;
end
subplot(2,1,2)
histfit(Stat)
[H,P,~] = swtest(Stat,0.05,0);
disp(['SW test fails to reject normality, P = ',num2str(P)])


figure('Position',[100,100,600,600])
histfit(Stat)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 1],'EdgeColor','k')
TrimFig(0.8);
SetFontSize(14);
xlabel('$$\frac{N_{k+1} - N_k}{N_k}$$','Interpreter','latex','Fontsize',15)
ylabel('bin count','Interpreter','latex','Fontsize',15)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600  ../../Figures/NormTest.png

figure('Position',[100,100,600,600])
normplot(Stat)
TrimFig(0.8);
SetFontSize(15);
title('')
set(gcf,'PaperPositionMode','auto')
print -dpng -r600  ../../Figures/NormTest2.png


%Levene Test for homoscedasticity
Year1 = Stat(1:52 - length(find((indexremoved >= 1) & (indexremoved <= 52))))';
i = length(Year1);
Year2 = Stat(i+1 : i+52 - length(find((indexremoved >= 53) & (indexremoved <= 104))))';
i = i + length(Year2);
Year3 = Stat(i+1 : i+52 - length(find((indexremoved >= 105) & (indexremoved <= 156))))';
i = i + length(Year3);
Year4 = Stat(i+1 :i+52 - length(find((indexremoved >= 157) & (indexremoved <= 208))))';
i = i + length(Year4);
Year5 = Stat(i+1 : i+52 - length(find((indexremoved >= 209) & (indexremoved <= 260))))';
i = i + length(Year5);
Year6 = Stat(i+1 : length(Stat) - length(find((indexremoved >= 261) & (indexremoved <= length(Stat)))))';

%ALL YEARS TOGETHER
Year1 = [Year1 1*ones(length(Year1),1)];
Year2 = [Year2 2*ones(length(Year2),1)];
Year3 = [Year3 3*ones(length(Year3),1)];
Year4 = [Year4 4*ones(length(Year4),1)];
Year5 = [Year5 5*ones(length(Year5),1)];
Year6 = [Year6 6*ones(length(Year6),1)];
X = [Year1;Year2;Year3;Year4;Year5;Year6];
Levenetest(X,0.05) % Rejects

%ONLY LAST 4 YEARS
Year3(:,2) = 1*ones(length(Year3),1);
Year4(:,2) = 2*ones(length(Year4),1);
Year5(:,2) = 3*ones(length(Year5),1);
Year6(:,2) = 4*ones(length(Year6),1);
X = [Year3;Year4;Year5;Year6];
Levenetest(X,0.05) % Fails to reject

%MLE of variance (see any text on stochastic volatility models)
MyData = Weeknum(104:end);
N = length(MyData);
Sigma2MLE = sum(log(MyData(2:end))-log(MyData(1:end-1)))/N - (log(MyData(end) - log(MyData(1))))^2/(N^3);
SigmaMLE = sqrt(Sigma2MLE)
