% Filename: ANSORepStudy.m
% Author: Andrew Zammit Mangion
% Date: November 2011
% Description: Studies the MC runs from RunMCwithcovariates.m
%
% Requires: Statistical Toolbox
%           histnorm.m (available from MATLAB Exchange)

clear all
close all

%ANSOCounts
load('GrowthtrialN2000WholePrecWithDistCovariates')
ANSOCount2009 = [43;65;48;88;292;137;10;100;23;85;239;52;227;1318;22;83;149;115;177;414;295;196;188;180;461;478;162;379;621;259;970;135];
ANSOCount2010 = [35;140;74;182;355;293;21;222;81;64;356;83;254;1457;4;84;195;126;146;511;504;356;263;491;1540;906;256;897;1387;353;1162;108];
Count2009(Count2009 == 0) = 1;
Count2010(Count2010 == 0) = 1;
ModelCount2009 = Count2009;
ModelCount2010 = Count2010;

% Simple linear predictor with no offset
Prediction = ModelCount2010./ModelCount2009.*repmat(ANSOCount2009,1,N);

figure('Position',[100,100,1200,500]); hold on;
for i = 1:length(Provname)
    Provname(i).strnum = strcat(num2str(i),'. ',Provname(i).str);
end
boxplot(log(Prediction'),{Provname(:).strnum},'labelorientation','inline'); hold on;
set(findobj(gca,'Type','text'),'FontSize',7)

% We are analyzing the log of the counts to stabilize the variances
plot(log(ANSOCount2010),'go','Linewidth',2)
Prediction(Prediction == 0) = 1;
Prediction(Prediction == Inf) = 1;
meanlogPred = mean(log(Prediction'))';
plot(meanlogPred,'ko','Linewidth',2)
ylabel('log AOG activity in 2010')
% set(gca,'YScale','log')
[r,p] = corrcoef(ANSOCount2010,median(Prediction')')
[r2,p2] = corrcoef(median(log(Prediction'))',log(ANSOCount2010))
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Predictions.png

% Point prediction comparison
figure('Position',[100 100 600 400]); hold on;
plot([1;8],[1;8],'r','Linewidth',2)
scatter(log(ANSOCount2010),mean(log(Prediction)')','k.','Linewidth',2)
xlabel('log ANSO count 2010')
ylabel('log model prediction 2010')
print -dpng -r300 Linearprediction.png

% Label points
figure('Position',[100 100 800 400]); hold on;
plot([1;8],[1;8],'r','Linewidth',2)
hold on
for i = 1:length(Provname)
    plot(log(ANSOCount2010(i)),median(log(Prediction(i,:)))','ko','Linewidth',2)
end
% Now get the points in the order in which they were plotted
ln = flipud(get(gca,'children'));  % ln(1) was plotted first.
% Now label them.
set(ln,'markers',16,'markerfa','w')
for ii = 2:length(ln)
    if (ii ~= 4) && (ii ~= 21) && (ii ~= 9) && (ii ~= 19) && ...
            (ii ~= 6) && (ii ~= 15) && (ii ~= 7)
        
        D = get(ln(ii),{'xdata','ydata'});
        S = sprintf('%i',ii-1);  % Second marker is first province
        T = text(D{1}-.04 * length(S),D{2},S) ;
        set(T,'fontsize',8,'color','k')
    end
end
text(3,5,'3');text(4,6,'18'); text(6,3.5,'8'); text(7,5,'20');
text(7,3.5,'6');text(5,7,'5'); text(8,6.5,'14');
xlabel('log ANSO count 2010')
ylabel('log model prediction 2010')
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Linearprediction_need_modif.png


% Point prediction with confidence intervals
figure('Position',[100 100 600 400]); hold on;
plot([1;8],[1;8],'r','Linewidth',2)
errorbar(log(ANSOCount2010),log(median(Prediction')'),...
    log(median(Prediction')')-log(quantile(Prediction',0.25)'),...
    log(quantile(Prediction',0.75)') - log(median(Prediction')'),'ko','Linewidth',2)
xlabel('log ANSO count 2010')
ylabel('log model prediction 2010')
print -dpng -r300 Linearprediction_withunc.png

% Histograms of log values
figure('units','normalized','outerposition',[0 0 1 1])
Pred_growth =  ((log(Prediction) - repmat(log(ANSOCount2009),1,N)))./(repmat(log(ANSOCount2009),1,N));
Act_growth = (log(ANSOCount2010) - log(ANSOCount2009))./log(ANSOCount2009);
Perc_change = [Act_growth,quantile(Pred_growth',[0.1,0.5,0.9])'];
for i = 1:length(Provname)
    subplot(6,6,i)
    histnorm(Pred_growth(i,:)*100,21); hold on
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k')
    %    histfit(Pred_growth(i,:)); hold on
    stem(Act_growth(i)*100,0.05);
    stem(median(Pred_growth(i,:))*100,0.05,'r');
    title(Provname(i).str,'Fontsize',8)
    axis([-70 70 0 0.05])
    set(gca,'Fontsize',8)
end
text(170,-0.01,'% growth in log activity count')
text(-320,0.18,'Normalised amplitude','rotation',90)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Growth_histogram.png

xlswrite('Perc_change.xls',{Provname.str}','Names')
xlswrite('Perc_change.xls',Perc_change,'Data')

% Kolmogorov-Smirnoff-like test
figure('Position',[100,100,700,400])
quantile_test = 0:0.0001:1; %Any smaller step size is noisy
class =zeros(size(quantile_test));

for i = 1:length(Provname)
    prov_quantiles = quantile(Pred_growth(i,:),quantile_test);
    closest_quant_ind = find(abs(prov_quantiles - Act_growth(i)) == min(abs(prov_quantiles - Act_growth(i))));
    closest_quant_ind =  closest_quant_ind(round(length(closest_quant_ind)/2));
    class(abs(closest_quant_ind(1) - find(quantile_test == 0.5))*2) = class(abs(closest_quant_ind(1) - find(quantile_test == 0.5))) + 1;
    text(quantile_test(closest_quant_ind),quantile_test(closest_quant_ind),Provname(i).str,'Fontsize',5)
end

plot(quantile_test,cumsum(class)/length(Provname),'k','Linewidth',2); hold on
plot([0 1],[0 1],'k--','Linewidth',2)
xlabel('quantiles'); ylabel('cumulative score'); box off
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Cumulative.png

% Save data
xlswrite('Perc_change.xls',{Provname.str}','Names')
xlswrite('Perc_change.xls',Perc_change,'Data')

