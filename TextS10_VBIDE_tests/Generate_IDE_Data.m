% Filename: Generate_IDE_Data.m
% Author: Andrew Zammit Mangion
% Date: April 2012
%
% Description: Generates PP data from underlying spatiotemporal AR1 process
%
% Requirements: STATISTICAL TOOLBOX
%               ../Common_functions (for PACF and PCCF computation)
%               poisson2d.m (generate 2D Homogeneous PP)



function Generate_IDE_Data(filename)

addpath('../Common_functions')
Anoise_low = 0.5;
Anoise_high = 2;
lnoise = 3;
bias = -2;
rho = 0.1;

%Setup time
dt = 1;
t = 0:dt:200; 
N = length(t);

%Setup space
J = 25; 
s = linspace(0,18,J);
ds = (max(s)-min(s))/(J-1);
[s1,s2] = meshgrid(s,s);

%Define Noise statistics
gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
Avec = repmat(Anoise_low,1,length(s));
Avariance_map = repmat(Avec,length(s),1);      %Atrue = 0.6

Avariance_map = Avariance_map + Anoise_high*gaussf(s1,s2,5,2.5,6,12);
Avariance_map = Avariance_map + Anoise_high*gaussf(s1,s2,10,13,6,12);

NoiseKernel.mu1 = (s(end)-s(1))/2 + s(1);
NoiseKernel.mu2 = (s(end)-s(1))/2 + s(1);
NoiseKernel.K = gaussf(s1,s2,NoiseKernel.mu1,NoiseKernel.mu2,lnoise,lnoise);

% Define variance field and growth field
s1vec = reshape(s1,[],1);
s2vec = reshape(s2,[],1);

% Vectorize displacements
[S1a,S1b] = meshgrid(s1vec);
[S2a,S2b] = meshgrid(s2vec);
S1a = S1a(:); S1b = S1b(:);
S2a = S2a(:); S2b = S2b(:);

%Generate variance field
Variance_cont =  @(s1,s2) interp2(Avariance_map,s1./ds+1,s2./ds+1);

%Generate growth field
K_growth = @(x1,x2) 0.2*exp(-((x1(:,1) - x2(:,1)).^2 + (x1(:,2) - x2(:,2)).^2)./lnoise);
Kvec2 = K_growth([S1a,S2a],[S1b,S2b]);
K2 = reshape(Kvec2,J^2,J^2);
K2 = triu(K2) + triu(K2)' - diag(diag(K2)); %Make positive definite and symmetric
z0 = zeros(J^2,1);
z = z0 + mvnrnd(z0,K2)';
Growth = reshape(z,J,J)*4;

% Generate noise sequences
Kf = @(x1,x2) Variance_cont(x1(:,1),x1(:,2)).*Variance_cont(x2(:,1),x2(:,2))...
          .*exp(-((x1(:,1) - x2(:,1)).^2 + (x1(:,2) - x2(:,2)).^2)./lnoise);
Kvec = Kf([S1a,S2a],[S1b,S2b]);
K = reshape(Kvec,J^2,J^2);
K = triu(K) + triu(K)' - diag(diag(K)); %Assure positive definite and symmetric
z0 = zeros(J^2,1);

noise = zeros(J,J,N);
for i = 1:N
    z = z0 + mvnrnd(z0,K)';
    noise(:,:,i) = reshape(z,J,J);
    i
end

%Evolve field
field = zeros(J,J,N);
lambda(:,:,1) = exp(bias + field(:,:,1));
for i = 2:N
    field(:,:,i) = rho*field(:,:,i-1)+Growth+ noise(:,:,i);
    lambda(:,:,i) = exp(bias + field(:,:,i));
end
z = @(i,s1,s2) interp2(field(:,:,i),s1./ds+1,s2./ds+1);

clear spikeaccept
%Generate events
spikeaccept.Coords = 0;
spikeaccept = repmat(spikeaccept,N,1);
for i = 1:N
    i
    area = (s(end)-s(1))^2;
    lambdastar = max(max(lambda(:,:,i)));
    spikecoords = poisson2d(lambdastar,area);
    spikeaccept(i).Coords = zeros(size(spikecoords));
    
    %Thin events
    k=1;
    r = rand(length(spikecoords),1);
    Uatpoint= z(i,spikecoords(:,1),spikecoords(:,2));
    NormIntatpoint =  exp(bias + Uatpoint)/lambdastar;
    for j = 1:length(spikecoords)
        if r(j) < NormIntatpoint(j) %Restrict to 0->1. Not explicitly done in Murray because of Sigmoid
            spikeaccept(i).Coords(k,:) = spikecoords(j,:);
            k = k+1;
        end
    end
    spikeaccept(i).Coords(k:end,:) = [];
end

save(filename)

