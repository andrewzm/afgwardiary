%Spatiotemporal nonparametric intensity estimator

clear all
close all

load('Results_with_Cov_estimation')
j  = 1;

for sigmas1 = [0.1 0.2 0.3 0.5 0.8]
        for sigmat = [1 1.5 2 2.5];
            sigmas2 = sigmas1;
            GRID = ST_Kernel_estimation(spikeaccept,sigmas1,sigmas2,sigmat,s1(1,:),s2(:,1),'T');
            save(strcat('Kernel_result',num2str(j)),'GRID','sigmas1','sigmas2','sigmat')
            j =  j+1;
        end
end
