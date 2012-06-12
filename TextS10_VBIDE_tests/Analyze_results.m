%--------------------------------------------------------------------------
% Filename: Analyze_results.m
% Authors: Andrew Zammit Mangion
% Date: March 19 2012
% 
% Details: Analyses results from VBIDE.m as in Section S10 in Supp material
%
% Requires: Results_with_cov_estimation.mat
%           Results_without_cov_estimation
%           Kernel_result1.mat - Kernel_result29.mat
%--------------------------------------------------------------------------



% Analyze estimation WITH and WITHOUT covariance
load('Results_with_cov_estimation')
Estinfo_with_cov = Estinfo;
Thetainfo_with_cov = Thetainfo;
load('Results_without_cov_estimation')
Estinfo_without_cov = Estinfo;
Thetainfo_without_cov = Thetainfo;

c = flipud(hot);

% Plot composite graph
figure('Position',[0,0,1100,700])
subplot(2,3,1)
for i = 1:50
    scatter(spikeaccept(i).Coords(:,1),spikeaccept(i).Coords(:,2),'k.')
    hold on
end
axis('tight')
xlabel('s_1'); ylabel('s_2');

subplot(2,3,4)
for i =1:Basis.nx
    %     surf(Basis.phi(:,:,i)); shading interp; hold on
    contour(s1,s2,Basis.phi(:,:,i),[1 exp(-1^2/2)],'r'); hold on
    plot(Basis.mu1(i),Basis.mu2(i),'kx')
end
view(2)
xlabel('s_1'); ylabel('s_2'); axis('tight')

s_sub = linspace(0,18,25);
[s1_sub,s2_sub] = meshgrid(s_sub);

subplot(2,3,2)
surf(s1_sub,s2_sub,Growth); caxis([-4 6]); shading interp; view(2);
colorbar('East'); axis('tight'); colormap(c)
xlabel('s_1'); ylabel('s_2');

subplot(2,3,3)
surf(s1,s2,Avariance_map.^2); caxis([0 6]); shading interp; view(2);
colorbar('East'); axis('tight'); colormap(c)
xlabel('s_1'); ylabel('s_2');

subplot(2,3,5)
Me_Growth = 0;
for i = 1:Basis.nx
    Me_Growth = Me_Growth + Thetainfo_with_cov(m).thetaest(i)*Basis.phi(:,:,i);
end;
surf(s1,s2,Me_Growth); caxis([-4 6]); shading interp; view(2);
colorbar('East'); axis('tight'); colormap(c)
xlabel('s_1'); ylabel('s_2');

subplot(2,3,6)
Cov = inv(Thetainfo_with_cov(m).Meanprecmat);
Me_Cov = 0;
for i = 1:Basis.nx
    Me_Cov = Me_Cov + Cov(i,i)*Basis.phi(:,:,i);
end;
colorbar; surf(s1,s2,Me_Cov); caxis([0 6]);  shading interp; view(2);
colorbar('East'); axis('tight'); colormap(c)
xlabel('s_1'); ylabel('s_2');

set(gcf,'PaperPositionMode','auto')
print -dpng -r600  CompFigure_Results.png

[FieldError1,TotfieldError1] = FindMSE_from_weights(s,[Estinfo_with_cov.xestRTS],field,Basis,'T');
[FieldError2,TotfieldError2] = FindMSE_from_weights(s,[Estinfo_without_cov.xestRTS],field,Basis,'T');
disp(['MSE_lambda with cov est = ',num2str(mean(FieldError1(2:end)))])
disp(['MSE_lambda without cov est = ',num2str(mean(FieldError2(2:end)))])

[GrowthError1,TotgrowthError1] = FindMSE_from_weights(s,Thetainfo_with_cov(m).thetaest,Growth,Basis,'T');
[GrowthError2,TotgrowthError2] = FindMSE_from_weights(s,Thetainfo_without_cov(m).thetaest,Growth,Basis,'T');
disp(['MSE_growth with cov est = ',num2str(GrowthError1)])
disp(['MSE_growth without cov est = ',num2str(GrowthError2)])

[c_cov,q,qMSE_cov,qbias_cov] = FindQuantiles_from_weights(s,Estinfo_with_cov(2:end-1),field,Basis);
[c_no_cov,q,qMSE_no_cov,qbias_no_cov] = FindQuantiles_from_weights(s,Estinfo_without_cov(2:end-1),field,Basis);
disp(['mean(bias_q) with cov est = ',num2str(mean(mean(qbias_cov)))])
disp(['mean(bias_q) without cov est = ',num2str(mean(mean(qbias_no_cov)))])


plot(q,mean(c_no_cov'),'k-d')
hold on; plot(q,mean(c_cov'),'k-s')
plot([0 1],[0 1],'r')

for i = 1:29
    load(strcat('Kernel_result',num2str(i)));
    FieldError3(:,i) = FindMSE(s,GRID,field);
    i
end
disp(['Kernel estimation MSE = ',num2str(min(mean(FieldError3)))]);
