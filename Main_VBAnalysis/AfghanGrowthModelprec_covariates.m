%--------------------------------------------------------------------------
% Filename: AfghanGrowthModelprec_covariates.m
% Author: Andrew Zammit Mangion
% Date: November 2011
% Description: An offline VBEM algorithm for the AWD
%
% Requires: ProcessedWeekData.mat
%           ../Shapefiles
%           ../TextS6_BasisSelection/AFGBasis (selected basis)
%           netlab (see http://www.fizyka.umk.pl/netlab/)
%           inpolygons.m (available from Matlab exchange)
%           posdist.m  (available from Matlab exchange)
%           multiprod.m  (available from Matlab exchange)
%
% Generates: ../AfghanModelHighResWholePrec_Covariates.mat
%--------------------------------------------------------------------------

function AfghanGrowthModelprec_covariates()

% load the data, readily organised into weekly bins
load('ProcessedWeekData')

Estmu = 'n';            % No need to estimate "background" since this is implied from basis selection
Estprec = 'y';          % Estimate precision matrix (inv(SigmaW))
Estcovariates = 'y';    % Estimate covariates

%-----------------
% Setup constants
%-----------------
tic

%Setup time
dt = 1 ;
N = 313;
t = 0:dt:dt*(N-1);

%Setup space
J = 101;
s = linspace(0,36,J);
ds = (max(s)-min(s))/(J-1);
l=5;         %length of basis functions
smax = s(end)-l; %These are spatial limits on basis representation
smin = l;
[s1,s2] = meshgrid(s,s);

% Note: For analysis Afghanistan is roughly warped into a tightly fitting
% square with left-hand corner on the origin and length 36. This is set
% through the scale/shift parameters. This affine transformation is not
% necessary.

% Find population density
Afg_pop = shaperead('../Shapefiles/afg_population_CSO_2011_2012','UseGeoCoords',true);
Pop_density = zeros(size(s1));
S1 = s1./scale(1) + shift(1);
S2 = s2./scale(2) + shift(2);
for i = 1:length(Afg_pop)
    X = inpolygons(S1(:),S2(:),Afg_pop(i).Lon',Afg_pop(i).Lat');
    Pop_density(X == 1) = Afg_pop(i).T_P_BS/ Afg_pop(i).AREA_SQKM;
end

% Find distance to city
Afg_cities = shaperead('../Shapefiles/07_03_settlements','UseGeoCoords',true);
Filtered_cities = 0; j =1;
for i = 1:length(Afg_cities)
    if Afg_cities(i).TYPE == 2 ||  Afg_cities(i).TYPE == 3
        Filtered_cities(j,1) = Afg_cities(i).LON;
        Filtered_cities(j,2) = Afg_cities(i).LAT;
        j = j+1;
    end
end
Dist_to_city = zeros(size(S1));
for i = 1:size(S1,1)
    for j = 1:size(S1,2)
        Dist_to_city(i,j) = min(posdist(S1(1,j),S2(i,1),Filtered_cities(:,1),Filtered_cities(:,2)));
    end
end


%Arrange in struct for par. passing
Constants.dt = dt;
Constants.ds = ds;
Constants.N = N;
Constants.t = t;
Constants.s = s;
Constants.s1 = s1;
Constants.s2 = s2;
Constants.smax = smax;
Constants.smin = smin;
Constants.J = J;
Constants.Pop_density = Pop_density;
Constants.Dist_to_city = Dist_to_city;


%-------------------------------------------------------------------
% Field setup
%-------------------------------------------------------------------

%Field noise
sigmaW = 0.2; % Rough estimate from temporal analysis, slightly de-tuned
sigmaR = 2;   % Observation noise, only used for linear regression on first time point

%Basis decomposition
load('../TextS6_BasisSelection/AFGBasis')       % Load chosen set of basis
Basis.inner = gaussinner(s1,s2,Basis.phi);      %#ok<NODEF> % Find the Gram matrix
Basis.Basisvec = zeros(J^2,Basis.nx);           % Vectorise
for i = 1:Basis.nx
    Basis.Basisvec(:,i) = reshape(Basis.phi(:,:,i),[],1);
end
Constants.d = Basis.nx;

%Setting up matrices and state-space model.
FieldMatrices.PSIx = Basis.inner;
FieldMatrices.W = dt*eye(Basis.nx);
FieldMatrices.R = sigmaR^2*eye(Basis.nx);

% --------------------
% Initialise variables
% --------------------
x = 3*ones(Basis.nx,N);          % Hidden states
numiters = 50;

% -------------------------------------------------------
% Initialise field estimate
% -------------------------------------------------------
Estinfo(1).xestpost = zeros(Basis.nx,1);        % Filtered estimate
Estinfo(1).xestprior = zeros(Basis.nx,1);       % Prior (for starting the scg) estimate
Estinfo(1).xestRTS = zeros(Basis.nx,1);         % Posterior (smoothed) estimate
Estinfo(1).PKalman = zeros(Basis.nx,Basis.nx);  % Filtered covariance
Estinfo(1).PRTS = zeros(Basis.nx,Basis.nx);     % Smoothed covariance
Estinfo = repmat(Estinfo(1),1,N);

% -------------------------------------------------------
% Initialise parameter estimate
% -------------------------------------------------------
Thetainfo(1).thetaest = zeros(Constants.d,1);
Thetainfo(1).vartheta = 1000*eye(Constants.d);
Thetainfo(1).muest = -3.5;
Thetainfo(1).varmu = 0;
Thetainfo(1).precisionest = repmat(1/(sigmaW^2),Basis.nx,1);  % Diagonal elements
Thetainfo(1).varprecision = repmat(1,Basis.nx,1);   
Thetainfo(1).Meanprecmat = 1/(sigmaW^2)*eye(Basis.nx);  % Whole matrix
Thetainfo(1).best = [-0.01,-0.01]; % Covariate weights
Thetainfo(1).varb = [0 0];        

Thetainfo = repmat(Thetainfo(1),1,numiters);

% -------------------------------------------------------
% Point process parameters
% -------------------------------------------------------
beta = 1; % Dummy variable if set to 1
Constants.beta = beta;


%-----------------------------------------------
% Inference
%-----------------------------------------------

%Estimate initial condition
EstIntensity = DiggleContSpace(spikeAllWeekData(1).Coords,Constants); % Nonparametric estimate of init condition
y = Regress(EstIntensity,Basis,Thetainfo(1).muest ,beta);             % Regress onto chosen basis
% Now we assume that the regression is a prior estimate for a Kalman filter
% to do a one-step update
Estinfo(1) = KalmanFilter(y,Estinfo(1),Constants,FieldMatrices,Basis);

% Start VBEM (Maximum of 200 VBEM iterations)
for m = 2:200
    % VBE step
    Estinfo = VBEMSmoother_full_nonlinear(spikeAllWeekData,Constants,Estinfo,FieldMatrices,Basis,Thetainfo(m-1));
    % VBM step
    Thetainfo(m) = VBM(Estinfo,Thetainfo(m-1),FieldMatrices,Constants,Basis,spikeAllWeekData,Estmu,Estprec,Estcovariates);
    % Save temporary results for debugging
    save('LatestResults')
    
    %Breaking Condition
    if norm(Thetainfo(m).thetaest - Thetainfo(m-1).thetaest) < 0.005 ...
            && (abs(Thetainfo(m).muest - Thetainfo(m-1).muest) < 0.01) ...
            && max(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) < 1.01 ...
            && min(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) > 0.99 ...
            && abs(Thetainfo(m).best(1) -Thetainfo(m-1).best(1)) < 0.005 ...
            && abs(Thetainfo(m).best(2) -Thetainfo(m-1).best(2)) < 0.005;
        break
    end
end

% Save results
save('../AfghanModelHighResWholePrec_Covariates')
toc
%----------------------------------------------------------



function [Estinfo] = VBEMSmoother_full_nonlinear(spikes,Constants,Initinfo,FieldMatrices,Basis,Thetainfo)

% -------------------------------------------------------------
% Function VBEMFilter_full_nonlinear
% inputs:
% outputs:
% -------------------------------------------------------------

% Initialize
SigmaW = inv(Thetainfo.Meanprecmat);
beta = Constants.beta;
b = Thetainfo.best;
popd_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.Pop_density,s1,s2);
distc_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.Dist_to_city,s1,s2);

xestprior = zeros(Basis.nx,Constants.N);
xestpost = zeros(Basis.nx,Constants.N);
xbeta = zeros(Basis.nx,Constants.N);
xestRTS = zeros(Basis.nx,Constants.N);

Sigmaprior = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmapost = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmabeta = zeros(Basis.nx,Basis.nx,Constants.N);
SigmaRTS = zeros(Basis.nx,Basis.nx,Constants.N);

Sigmaprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
Sigmapost(:,:,1) = 30*eye(Basis.nx);
xestpost(:,1) = Initinfo(1).xestpost;

mu = Thetainfo.muest;
theta = Thetainfo.thetaest;
Qinv = inv(SigmaW);
dt = Constants.dt;
Estinfo = Initinfo;

% Options
options = foptions;
options(14) = 2000;
% options(9) = 1;
options(2) = 0.1;
options(3) = 0.1;

% Forward pass
for i = 2:Constants.N
    xestprior(:,i) = xestpost(:,i-1) + Constants.dt*theta;
    Sigmaprior(:,:,i) = Sigmapost(:,:,i-1) + SigmaW;
    
    Sigmatilde = inv(inv(Sigmapost(:,:,i-1)) + Qinv);
    Sigmastar = inv(Qinv - Qinv*Sigmatilde*Qinv);
    mustar = Sigmastar*(Qinv*Sigmatilde*(inv(Sigmapost(:,:,i-1))*xestpost(:,i-1) - Qinv*theta*dt) + Qinv*theta*dt);
    Sigmastarinv = inv(Sigmastar);
    
    spikecoords = spikes(i).Coords;
    
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        popd_eval = popd_cont(spikecoords(:,1),spikecoords(:,2));
        distc_eval = distc_cont(spikecoords(:,1),spikecoords(:,2));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.sigma2(j),Basis.sigma2(j))';
        end
    else
        phieval = zeros(Basis.nx,size(spikecoords,1));
        popd_eval = 0; distc_eval = 0;
    end
    
    myint = @(xx)  SingleIntForward(xx,Constants,Basis,Thetainfo);
    myint2 = @(xx)  MultiIntForward(xx,Constants,Basis,Thetainfo);
    f = @(xx) -(sum(mu + b(1)*popd_eval + b(2)*distc_eval + beta*phieval'*xx') - myint(xx') - ((xx' - mustar)'/Sigmastar)*(xx' - mustar)./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmastar + mustar'/Sigmastar);
    %[temp,options,t1,t2,t3] = scg(f,xestprior(:,i)',options,gradf);
    [xestpost(:,i)] = scg(f,xestprior(:,i)',options,gradf)';
    myint3 = @(xx) MultiIntForward2(xx,Constants,Basis,Thetainfo);
    temp = inv(Sigmastar) + myint3(xestpost(:,i));
    Sigmapost(:,:,i) = inv(temp);
    [i max(xestpost(:,i)) min(xestpost(:,i))]
    Estinfo(i).xestpost = xestpost(:,i);
    Estinfo(i).xestprior = xestprior(:,i);
    Estinfo(i).PKalman = Sigmapost(:,:,i);
end


% Backward pass
Sigmabeta(:,:,i) = 9*eye(Basis.nx);
xbeta(:,end) = xestpost(:,end);
for i = Constants.N-1:-1:1
    spikecoords = spikes(i+1).Coords;
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        popd_eval = popd_cont(spikecoords(:,1),spikecoords(:,2));
        distc_eval = distc_cont(spikecoords(:,1),spikecoords(:,2));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.sigma2(j),Basis.sigma2(j))';
        end
    else
        phieval = zeros(Basis.nx,size(spikecoords,1));
        popd_eval = 0; distc_eval = 0;
    end
    
    f = @(xx) -(sum(mu + b(1)*popd_eval + b(2)*distc_eval + beta*phieval'*xx') - myint(xx') - ((xx' - xbeta(:,i+1))'/Sigmabeta(:,:,i+1))*(xx' - xbeta(:,i+1))./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmabeta(:,:,i+1) + xbeta(:,i+1)'/Sigmabeta(:,:,i+1));
    mudash = scg(f,xestpost(:,i+1)',options,gradf)';
    Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(mudash));
    
    Sigmatilde = inv(inv(Sigmadash) + Qinv);
    Sigmabeta(:,:,i) = inv(Qinv - Qinv*Sigmatilde*Qinv);
    xbeta(:,i) = Sigmabeta(:,:,i)*(-dt*Qinv*theta + Qinv*Sigmatilde*(inv(Sigmadash)*mudash + Qinv*dt*theta));
    
    SigmaRTS(:,:,i) = inv(inv(Sigmapost(:,:,i)) + inv(Sigmabeta(:,:,i)));
    xestRTS(:,i) = SigmaRTS(:,:,i)*(inv(Sigmapost(:,:,i))*xestpost(:,i) + inv(Sigmabeta(:,:,i))*xbeta(:,i));
    
    i
    
    Estinfo(i).xestRTS = xestRTS(:,i);
    Estinfo(i).PRTS = SigmaRTS(:,:,i);
    
    temp = Qinv + inv(Sigmabeta(:,:,i+1))+ myint3(xestRTS(:,i+1));
    Sigmatilde = inv(inv(Sigmapost(:,:,i)) + Qinv);
    Crosscov = Sigmatilde*(Qinv)*inv(temp -  Qinv*Sigmatilde*Qinv);
    
    Estinfo(i).W = xestRTS(:,i)*xestRTS(:,i)' + SigmaRTS(:,:,i);
    Estinfo(i).S = xestRTS(:,i)*xestRTS(:,i+1)' + Crosscov;    
end

%Update posterior for next time step
xestpost(:,1) = xestRTS(:,2) - theta*dt; %initial state to the next iteration
% Sigmapost(:,:,1) = SigmaW*inv(eye(Basis.nx)-A*A'); %initial state

function [Estinfo] = KalmanFilter(y,Initinfo,Constants,FieldMatrices,Basis)


% -------------------------------------------------------------
% Function KalmanFilter
% inputs: y (prior), and ST structures.
% outputs: posterior on initial field
% Description: One-step update of the initial condition
% -------------------------------------------------------------

% Initialize
SigmaR = FieldMatrices.R;
xestpost = zeros(Basis.nx,Constants.N);
sigma2estpost = zeros(Basis.nx,Basis.nx,Constants.N);
xestprior = zeros(Basis.nx,Constants.N);
sigma2estprior = zeros(Basis.nx,Basis.nx,Constants.N);

sigma2estprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
xestpost(:,1) = Initinfo(1).xestpost;

Estinfo = Initinfo;
% Assume point observation
C = eye(Basis.nx);

S = C*sigma2estprior(:,:,1)*C' + SigmaR;
K = sigma2estprior(:,:,1)*C'*inv(S); % Kalman gain
xestpost(:,1) = xestprior(:,1) + K*(y(:,1) - C*xestprior(:,1)); % Mean update
Estinfo(1).xestpost = xestpost(:,1); % Return initial condition

function [lambda] = DiggleContSpace(spikes,Constants)

% -------------------------------------------------------------
% Function DiggleContSpace
% inputs: coordinates of events, Constants structure
% outputs: Nonparametric estimate of spatial intensity
% Note: This function is overly complicated, it would be better to 
% use standard kernels and then digitise at the end.
%---------------------------------------------------------------

% Initialise
lambda = zeros(Constants.J,Constants.J);
frame = lambda;

% A simple (quick) cylindrical estimator on a digitised map
buff = 8;
rabs = buff*Constants.ds;
r = buff;

% Find pixel where event lands 
for i = 1:size(spikes,1)
    %Find x pixel
    [~,x] = find((spikes(i,1) - Constants.s).^2 == min((spikes(i,1) - Constants.s).^2));
    %Find y pixel
    [~,y] = find((spikes(i,2) - Constants.s).^2 == min((spikes(i,2) - Constants.s).^2));
    frame(y,x) = frame(y,x)+1;
end

% Add a digitised (normalised) cylinder for every event in pixel
for i = buff:Constants.J-buff
    for j = buff:Constants.J-buff
        lambda(j,i) = pixincircle(frame,i,j,r)/(pi*rabs^2)/Constants.dt;
    end
end

function y = pixincircle(A,xc,yc,radius)

% smooth out (xc,yc) based on the count map A using cylinders of radius "radius" pixels 
% (finds how many circles from events enclose this (xc, yc) )

y=1:size(A,1);
x=1:size(A,2);
[X Y]=meshgrid(x,y);

% This assume the circle falls *entirely* inside the image
R2 = (X-xc).^2+(Y-yc).^2;
[c1,c2] = find(R2 < radius^2);
c(1,:) = c1;
c(2,:) = c2;
c = round(c(:,2:end)); % pixels located in circle
Ac = A(sub2ind(size(A),c(1,:),c(2,:))); % extract value
y = sum(Ac);


function y = Regress(lambda,Basis,mu,beta)

% Regress an intensity map "lambda" onto the selected basis using mean
% square.

% Initialize
phivec = (reshape(Basis.phi,[],Basis.nx))';
y = zeros(Basis.nx,size(lambda,3));
lambda = lambda + 0.001;  % Some regularisation

% Now find the LS estimate of the initial field
for i = 1:size(lambda,3)
    lambdavec = reshape(lambda(:,:,i),[],1);
    y(:,i) = inv(phivec*phivec')*phivec*((log(lambdavec) - mu)/beta);
end

function [Thetainfo] = VBM(Estinfo,Thetainfo,FieldMatrices,Constants,Basis,spikes,Estmu,Estprec,Estcovariates)

SigmaW = diag(1./Thetainfo.precisionest)*FieldMatrices.W;
d = Constants.d;
dt = Constants.dt;
N = Constants.N;



%Estimate theta
v = zeros(d,1);
for i = 3:Constants.N-2
    v(:,1) = v(:,1) + dt*inv(SigmaW)*(Estinfo(i).xestRTS - Estinfo(i-1).xestRTS);
end
vartheta = SigmaW./dt^2/(Constants.N-4);
thetaest = vartheta*v;
Thetainfo.thetaest = thetaest;
Thetainfo.vartheta = vartheta;

%Estimate precision
alpha = 500;
beta = 20;

precest = repmat(1/(0.2^2),Basis.nx,1);
varprec = repmat(0,Basis.nx,1);
if Estprec == 'y'
    % Diagonal estimation
    for j = 1:Basis.nx
        k1 = 0;
        for i = 3:N-2
            k1 = k1+ Estinfo(i).W(j,j) + Estinfo(i-1).W(j,j) + dt^2*(vartheta(j,j) + thetaest(j)^2) ...
                - 2*Estinfo(i-1).S(j,j) + 2*dt*thetaest(j)*Estinfo(i-1).xestRTS(j) ...
                - 2*dt*thetaest(j)*Estinfo(i).xestRTS(j);
        end
        alpha0 = alpha + (Constants.N-4)/2;
        beta0 = beta + k1/2;
        precest(j) = alpha0/beta0;
        varprec(j) = alpha0/beta0^2;
    end
    Meanprecmat = diag(precest);
    
    %Whole Wishart distribution
    dofprior = 1000;
    Prevmatprior = 25/dofprior*eye(Basis.nx);
    Gamma = 0;
    for i = 3:N-2
        Gamma = Gamma + Estinfo(i).W + Estinfo(i-1).W + vartheta + thetaest*thetaest' ...
            - Estinfo(i-1).S - Estinfo(i-1).S' ...
            - Estinfo(i).xestRTS*thetaest' - thetaest* Estinfo(i).xestRTS'  ...
            + Estinfo(i-1).xestRTS*thetaest' + thetaest*Estinfo(i-1).xestRTS';
    end
    Precmatpost = inv(inv(Prevmatprior) + Gamma);
    Meanprecmat = (dofprior + N-4)*Precmatpost;
end

%Estimate mu (Convergence not assured due to identifiability)
if Estmu == 'n'
    Thetainfo(1).muest = -3.5;
    Thetainfo(1).varmu = 0;
else
    muprior = 0;
    Sigmamuprior = 10;
    totalspikenum = 0;
    totalvoidsum = 0;
    phitemp = zeros(Constants.J,Constants.J,Basis.nx);
    for i = 3:N-2
        Umean = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS,1,1,Basis.nx)),3);
        MySigma = Estinfo(i).PRTS; %CHANGE TO FULL COVARIANCE
        for j = 1:Basis.nx
            phitemp(:,:,j) = sum(multiprod(Basis.phi,reshape(MySigma(:,j),1,1,Basis.nx)),3);
        end
        phisigma = sum(phitemp.*Basis.phi,3);
        totalvoidsum = totalvoidsum + dt*Constants.ds^2*sum(sum(exp(beta*Umean + beta^2*phisigma./2)));
        totalspikenum = totalspikenum + size(spikes(i).Coords,1);
    end
    fmu = @(mu) -(-(mu - muprior)^2/(2*Sigmamuprior) + mu*totalspikenum - exp(mu)*totalvoidsum);
    gradfmu = @(mu) -(-mu/Sigmamuprior + totalspikenum - exp(mu)*totalvoidsum);
    options = foptions;
    % options(9) = 1;
    Thetainfo(1).muest = scg(fmu,muprior,options,gradfmu)';
    Thetainfo(1).varmu = 1/(1/Sigmamuprior + exp(Thetainfo(1).muest)*totalvoidsum);
end

%Estimate covariates
if Estcovariates == 'n'
    Thetainfo(1).b(1) = 0;
    Thetainfo(1).b(2) = 0;
    Thetainfo(1).varb(1) = 0;
    Thetainfo(1).varb(2) = 0;
else
    bprior(1) = 0;
    bprior(2) = 0;
    Sigmabprior(1) = 10;
    Sigmabprior(2) = 10;
    pop_den_sum = 0;
    distc_sum = 0;
    totalfieldsum = 0;
    popd_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.Pop_density,s1,s2);
    distc_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.Dist_to_city,s1,s2);
    Basis_vec_mat = reshape(Basis.phi,[],Basis.nx);
    phitemp = zeros(Constants.J,Constants.J,Basis.nx);
    for i = 3:N-2
        Umean = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS,1,1,Basis.nx)),3);
        MySigma = Estinfo(i).PRTS;
        phitemp_vec = Basis_vec_mat*MySigma;
        phitemp = reshape(phitemp_vec,Constants.J,Constants.J,Basis.nx);
        phisigma = sum(phitemp.*Basis.phi,3);
        totalfieldsum = totalfieldsum + dt*(exp(Constants.beta*Umean + Constants.beta^2*phisigma./2));
        pop_den_sum = pop_den_sum + sum(popd_cont(spikes(i).Coords(:,1),spikes(i).Coords(:,2)));
        distc_sum = distc_sum + sum(distc_cont(spikes(i).Coords(:,1),spikes(i).Coords(:,2)));
    end
    b1int_1 = @(b1) sum(sum(exp(Thetainfo.muest +Thetainfo.varmu/2+ b1*Constants.Pop_density + Thetainfo.best(2)*Constants.Dist_to_city + Thetainfo.varb(2)*Constants.Dist_to_city.^2./2).*totalfieldsum))*Constants.ds^2;
    b1int_2 = @(b1) sum(sum(Constants.Pop_density.*exp(Thetainfo.muest + b1*Constants.Pop_density + Thetainfo.best(2)*Constants.Dist_to_city +Thetainfo.varb(2)*Constants.Dist_to_city.^2./2).*totalfieldsum))*Constants.ds^2;
    b1int_3 = @(b1) sum(sum((Constants.Pop_density.^2).*exp(Thetainfo.muest + b1*Constants.Pop_density + Thetainfo.best(2)*Constants.Dist_to_city + Thetainfo.varb(2)*Constants.Dist_to_city.^2./2).*totalfieldsum))*Constants.ds^2;
    
    fb1 = @(b1) -(-(b1 - bprior(1))^2/(2*Sigmabprior(1)) + b1*pop_den_sum - b1int_1(b1));
    gradfb1 = @(b1) -(-(b1 - bprior(1))/Sigmabprior(1) + pop_den_sum - b1int_2(b1));
    options = foptions;
    options(9) = 1;
    Thetainfo(1).best(1) = scg(fb1,bprior(1),options,gradfb1)';
    Thetainfo(1).varb(1) = 1/(1/Sigmabprior(1) + b1int_3(Thetainfo(1).best(1)));
    
    b2int_1 = @(b2) sum(sum(exp(Thetainfo.muest + Thetainfo.best(1)*Constants.Pop_density + Thetainfo.varb(1)*Constants.Pop_density.^2/2 + b2*Constants.Dist_to_city).*totalfieldsum))*Constants.ds^2;
    b2int_2 = @(b2) sum(sum(Constants.Dist_to_city.*exp(Thetainfo.muest + Thetainfo.best(1)*Constants.Pop_density + Thetainfo.varb(1)*Constants.Pop_density.^2/2 + b2*Constants.Dist_to_city).*totalfieldsum))*Constants.ds^2;
    b2int_3 = @(b2) sum(sum((Constants.Dist_to_city.^2).*exp(Thetainfo.muest + Thetainfo.best(1)*Constants.Pop_density + Thetainfo.varb(1)*Constants.Pop_density.^2/2 + b2*Constants.Dist_to_city).*totalfieldsum))*Constants.ds^2;
    
    fb2 = @(b2) -(-(b2 - bprior(2))^2/(2*Sigmabprior(2)) + b2*distc_sum - b2int_1(b2));
    gradfb2 = @(b2) -(-(b2 - bprior(2))/Sigmabprior(2) + distc_sum - b2int_2(b2));
    Thetainfo(1).best(2) = scg(fb2,bprior(2),options,gradfb2)';
    Thetainfo(1).varb(2) = 1/(1/Sigmabprior(2) + b2int_3(Thetainfo(1).best(2)));
end


Thetainfo(1).precisionest = precest;
Thetainfo(1).varprecision = varprec;
Thetainfo(1).Meanprecmat = Meanprecmat;

function [phi] = LocalisedKernelPhi_Cont(s1,s2,mu1,mu2,sigma21,sigma22)

% Evaluate CGRBF with centre (mu1,mu2) with stds (sigma21,sigma22) at
% points (s1,s2)

phi = zeros(length(s1),1);
beta1 = sqrt(pi/sigma21);
beta2 = sqrt(pi/sigma22);
l1 = 2*pi/beta1;
l2 = 2*pi/beta2;
slow1 = (mu1 - l1);  %Find limit kernel
shigh1 = (mu1 + l1);  %Find limit kernel
slow2 = (mu2 - l2);  %Find limit kernel
shigh2 = (mu2 + l2);  %Find limit kernel

Delta1 = beta1*abs(s1 - mu1);
Delta2 = beta2*abs(s2 - mu2);
phi(:,1) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
phi = phi.*~(s1 < slow1 | s1 > shigh1 | s2 < slow2 | s2 > shigh2);


function y = gaussinner(s1,s2,phi)
% Find the inner product on s1,s2 of a 3-D array of basis functions
nx = size(phi,3);
for i = 1:nx
    for j = 1:nx
        y(i,j) = trapz(s2(:,1),(trapz(s1(1,:),phi(:,:,i).*phi(:,:,j),2)));
    end
end

function y = SingleIntForward(x,Constants,Basis,Thetainfo)
% An integral required in optimisation of VB-Laplace
beta = Constants.beta;
mu = Thetainfo.muest;
b = Thetainfo.best;
varb = Thetainfo.varb;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y = Constants.dt*Constants.ds^2*exp(mu + varmu/2)*sum(exp(b(1)*Constants.Pop_density(:) + varb(1)*Constants.Pop_density(:).^2./2  ...
    + b(2)*Constants.Dist_to_city(:) + varb(2)*Constants.Dist_to_city(:).^2./2+ beta*Uest));

function y = MultiIntForward2(x,Constants,Basis,Thetainfo)
% An integral required in optimisation of VB-Laplace
beta = Constants.beta;
mu = Thetainfo.muest;
b = Thetainfo.best;
varb = Thetainfo.varb;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y  = Constants.dt*Constants.ds^2*exp(mu+varmu/2).*(beta^2*Basis.Basisvec'*(Basis.Basisvec.*repmat(exp(b(1)*Constants.Pop_density(:) + varb(1)*Constants.Pop_density(:).^2./2  ...
    + b(2)*Constants.Dist_to_city(:) + varb(2)*Constants.Dist_to_city(:).^2./2 + beta*Uest),1,Basis.nx)));

function y = MultiIntForward(x,Constants,Basis,Thetainfo)
% An integral required in optimisation of VB-Laplace
beta = Constants.beta;
mu = Thetainfo.muest;
b = Thetainfo.best;
varb = Thetainfo.varb;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y =  Constants.dt*Constants.ds^2*exp(mu + varmu/2).*(beta*Basis.Basisvec'*exp(beta*Uest + b(1)*Constants.Pop_density(:) +  varb(1)*Constants.Pop_density(:).^2./2 ...
    + b(2)*Constants.Dist_to_city(:) +  varb(2)*Constants.Dist_to_city(:).^2./2));


