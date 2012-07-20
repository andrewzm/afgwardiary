% Filename: AfghanCorrelationTestsPlot.m
% Author: Andrew Zammit Mangion
% Date: June 2011
%
% Requires: Image Processing Toolbox
%           CorrelationTests.mat
%           ../Common_functions
%           jbfill.m (Available from Matlab Exchange)
% Description: Computes the noise kernel and mixing kernel of the homogeneous
% IDE. Finally, from the PACF, the parameters for the basis (tau or/and sigma2b) are
% found.

clear all
close all

addpath('../Common_functions')
load('CorrelationTests')

%Setup space
J = 101;
s = linspace(0,36,J);
ds = (max(s)-min(s))/(J-1);
l=5;                % support length of basis functions
smax = s(end)-l;    % These are spatial limits on basis representation
smin = l;
[s1,s2] = meshgrid(s,s);


gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
mu1 = (s(end)-s(1))/2 + s(1);
mu2 = (s(end)-s(1))/2 + s(1);
savefigures = 'n';

%Parametric estimation of log(g_{k,k}) and log(g_{k,k+1})
[Anoptim,lnoptim] = MethodContrast(r(2:length(r)),exp(mean(log(gest(2:end,2:end)),2)),'RBF',r(2),r(length(r)));
ObservedCovMean = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
[Akoptim,lkoptim] = MethodContrast(r(3:end),exp(mean(log(gcrossest(3:end,2:end)),2)),'RBF',r(3),r(end));
ObservedCrossCovMean = Akoptim*gaussf(s1,s2,mu1,mu2,lkoptim,lkoptim);

%Nonparametric estimation of IDE Kernel using deconvolution with some
%regularization
Kernelest = real(deconvwnr(ObservedCrossCovMean,ObservedCovMean,0.0000001))/ds^2;

%Parametric reduction
Koptim = @(s1,s2,par) ds^2*sum(sum((par(1)*exp(-(s1-mu1).^2/par(2)-(s2-mu2).^2/par(2)) - Kernelest).^2));
parest = fminsearch(@(par) Koptim(s1,s2,par),[1,1]);
Anoptim = parest(1);
lnoptim = parest(2);
Kernelest = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
midind = ceil(J/2);


%Plot kernel
f = figure('Position',[100,100,700*0.75,500*0.75]);
plot(s-s(midind),diag(Kernelest),'k','Linewidth',2)
xlabel('\upsilon')
ylabel('k_I(\upsilon)')
SetFontSize(18);
TrimFig(0.9)
axis([-18 18 0 0.7])
set(gcf,'PaperPositionMode','auto')
if savefigures == 'y'
        print -dpng -r600  Kernel.png
end

%Plot ln g_k,k
f = figure('Position',[100,100,700*0.75,500*0.75]);
plot(s-s(midind),diag(ObservedCovMean),'k','Linewidth',2); hold on
diagCrossCovMean = diag(ObservedCrossCovMean);
plot(s(1:2:end)-s(midind),diagCrossCovMean(1:2:end),'kx','Linewidth',1.5)
xlabel('\upsilon')
ylabel('\cdot(\upsilon)')
SetFontSize(18);
TrimFig(0.9)
% leg = legend('$$\ln~\bar{g}_{k,k}~~~(\upsilon)$$','$$\ln~\bar{g}_{k,k+1}~~~~(\upsilon)$$')
% set(leg,'interpreter','latex','Fontsize',7)
set(gcf,'PaperPositionMode','auto')
axis([-15,15,0,2.8])
if savefigures == 'y'
        print -dpng -r600  AutoCov.png
end


f = figure('Position',[100,100,600,500]);
[Anoptim,lnoptim] = MethodContrast(r(2:end),exp(mean(log(gest(2:end,1:150)),2)),'RBF',r(2),r(end));
InitObservedCovMean = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
[Anoptim,lnoptim] = MethodContrast(r(2:end),exp(mean(log(gest(2:end,150:end)),2)),'RBF',r(2),r(end));
FinalObservedCovMean = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
% plot(s(1:midind-12),EstNoiseCov(midind:end-12,midind),'r','LineWidth',2)           %Estimated noise Kernel
plot(s(1:midind-12),ObservedCovMean(midind:end-12,midind)/max(max(ObservedCovMean)),'m','LineWidth',2); hold on           %Estimated noise Kernel
plot(s(1:midind-12),InitObservedCovMean(midind:end-12,midind)/max(max(InitObservedCovMean)),'g--','LineWidth',2)           %Estimated noise Kernel
plot(s(1:midind-12),FinalObservedCovMean(midind:end-12,midind)/max(max(FinalObservedCovMean)),'b--','LineWidth',2)           %Estimated noise Kernel
plot(s,exp(-(s).^2/(2*0.8106)),'k-.','LineWidth',2)           %Estimated noise Kernel
xlabel('\upsilon','FontSize',12)
axis([0 15 0 1])
SetFontSize(32);
leg = legend('$$\ln\bar{g}_{k,k}~(\upsilon)$$ (2004-2009)','$$\ln\bar{g}_{k,k}~(\upsilon)$$  (2004-2006)','$$\ln\bar{g}_{k,k}~(\upsilon)$$  (2007-2009)','$$\phi(\upsilon)$$')
set(leg,'Interpreter','latex','Fontsize',20)
set(leg,'Position',get(leg,'Position') - [0.02 0 0 0])
legend('boxoff')
[Anoptim,lnoptim] = MethodContrast(r(2:end),exp(mean(log(gest(2:end,1:end)),2) + 2*std(log(gest(2:end,1:end))')'),'RBF',r(2),r(end));
UpperLimit = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
[Anoptim,lnoptim] = MethodContrast(r(2:end),exp(mean(log(gest(2:end,1:end)),2) - 2*std(log(gest(2:end,1:end))')'),'RBF',r(2),r(end));
LowerLimit = Anoptim*gaussf(s1,s2,mu1,mu2,lnoptim,lnoptim);
jbfill(s(1:midind-12),UpperLimit(midind:end-12,midind)'/max(max(UpperLimit)),LowerLimit(midind:end-12,midind)'/max(max(LowerLimit)),'y');
TrimFig(0.9);
set(gca, 'Box', 'off' );
if savefigures == 'y'
        print -dpng -r600  ../../Figures/AutoCovAll.png
end


%The PACF is proportional to the autocorrelation function
%The Fourier transform of the autocorrelation function is the PSD
%sqrt of the PSD is the freq. response (Zammit_2012)
f = figure('Position',[100,100,600,500]);
Cutoff = 0.2;
[S] =0:J-1;
S = S./(J*ds);
Freqresponse =  sqrt(abs(fft(diag(ObservedCovMean))));
plot(S,Freqresponse./Freqresponse(1),'b','Linewidth',2); hold on
plot(S,exp(-(S).^2/(Cutoff^2)),'r','Linewidth',2); hold on
plot([Cutoff Cutoff], [0 1], 'k','Linewidth',2)
xlabel('\nu')
axis([0 0.6 0 1])
SetFontSize(32)
ylabel('$$|\mathcal{F}(\cdot)|$$','Interpreter','latex','Fontsize',32)
TrimFig(0.9)
set(gca, 'Box', 'off' );
if savefigures == 'y'
        print -dpng -r600  ../../Figures/FreqResponse.png
end

BasisVar = 1/(2*Cutoff^2*pi^2)
tau = sqrt(2*Cutoff^2*pi^3)
