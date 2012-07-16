%--------------------------------------------------------------------------
% Program VBEM for IDE
% Authors: Andrew Zammit Mangion
% Date: March 19 2012
% 
% Details: An offline VBEM algorithm for state (field) + parameter
% inference of a spatiotemporal state-space system with point process observations
%
% Input: filename (output filename)
%        Estprec ('y' or 'n'. If 'n' the identity matrix is used for the precision)
%--------------------------------------------------------------------------

function y = VBIDE(filename,Estprec)

load('Test_data')
 
Esttheta = 'y';
Estkernel = 'n';
Estrho = 'n';
Estmu = 'n';

%Refine space
J = 101;
s = linspace(s(1),s(end),J); %#ok<NODEF>
ds = (max(s)-min(s))/(J-1);
% [tau,cutoff]=FrequencyAnalysis(r,gest,s);         %length of basis functions 0.15 for Fine data
tau = 1.4175;
cutoff= 0.18;
smax = s(end);
smin = 0;
[s1,s2] = meshgrid(s,s);


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


%-------------------------------------------------------------------
% Field setup
%-------------------------------------------------------------------

%Field noise
sigmaW = 0.2;
sigmaR = 2;

spacing = 1/(2*cutoff*1.2);
mu1 = linspace(s(1),s(end),(s(end)-s(1))/spacing + 1);
mu2 = linspace(s(1),s(end),(s(end)-s(1))/spacing + 1);
Basis.nx = length(mu1)^2;

% Field basis functions
[C1] = meshgrid(mu1,mu2);
Basis.mu1 = reshape(C1,1,[]);
Basis.mu2 = reshape(C1',1,[]);
tau1(1:Basis.nx) = tau;
tau2(1:Basis.nx) = tau;
Basis.tau1 = tau1;
Basis.tau2 = tau2;
Basis.phi = LocalisedKernelPhi(s1,s2,Basis.mu1,Basis.mu2,tau1,tau2);
close all
for i =1:Basis.nx
    surf(Basis.phi(:,:,i)); shading interp; hold on
end

if Estkernel == 'y'
    % Kernel basis functions
    gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
    KernelBasis.d = 1;
    KernelBasis.phi(:,:,1) = gaussf(s1,s2,IDEKernel.mu1,IDEKernel.mu2,IDEKernel.sigma2,IDEKernel.sigma2);
end

% Required matrices
Basis.inner = gaussinner(s1,s2,Basis.phi);
Basis.Basisvec = zeros(J^2,Basis.nx);
for i = 1:Basis.nx
    Basis.Basisvec(:,i) = reshape(Basis.phi(:,:,i),[],1);
end
if Estkernel == 'y'
    KernelBasis.PHI = zeros(length(s),length(s),KernelBasis.d,Basis.nx);
    for j = 1:KernelBasis.d
        for i = 1:Basis.nx
            KernelBasis.PHI(:,:,j,i) = conv2(Basis.phi(:,:,i),KernelBasis.phi(:,:,j),'same')*ds^2;
        end
    end
end
%Noise mean
b(1,1:Basis.nx) = 0;
Constants.d = length(b);

%Setting up matrices and state-space model.
FieldMatrices.PSIx = Basis.inner;

% Function definitions
[S1x,S2x] = meshgrid(linspace(s(1),s(end),length(Growth)));
gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
Growth_interp = @(s1,s2) interp2(S1x,S2x,Growth,s1,s2);
Avariance_cont =  @(s1,s2) interp2(S1x,S2x,Avariance_map,s1,s2); %#ok<NODEF>

%Exact reduced noise covariance matrix
NoiseKernel.K = gaussf(s1,s2,NoiseKernel.mu1,NoiseKernel.mu2,lnoise,lnoise); %#ok<NODEF>
Avariance_map = Avariance_cont(s1,s2);
Qphi = zeros(J,J,Basis.nx);
for i = 1:Basis.nx
    Qphi(:,:,i) = conv2(Avariance_map.*Basis.phi(:,:,i),NoiseKernel.K,'same')*ds^2;
end
FieldMatrices.Qn = gaussinner2(s1,s2,repmat(Avariance_map,[1,1,Basis.nx]).*Basis.phi,Qphi);
W = inv(FieldMatrices.PSIx)*FieldMatrices.Qn*inv(FieldMatrices.PSIx);
W2 = triu(W) + triu(W)' - diag(diag(W)); %Ensure symmetry!! (numerical errors)
FieldMatrices.W = W2;

if Estkernel == 'y' 
    FieldMatrices.V = FindV(Constants.s1,Constants.s2,Basis,KernelBasis);
end

%Growth vector
Basis_vec_mat = reshape(Basis.phi,[],Basis.nx);
thetatrue = inv(Basis_vec_mat'*Basis_vec_mat)*Basis_vec_mat'*reshape(Growth_interp(s1,s2),[],1);


% --------------------
% Initialise variables
% --------------------
numiters = 50;

% -------------------------------------------------------
% Initialise field estimate
% -------------------------------------------------------
Estinfo(1).xestpost = zeros(Basis.nx,1);
Estinfo(1).xestprior = zeros(Basis.nx,1);
Estinfo(1).xestRTS = zeros(Basis.nx,1);
Estinfo(1).PKalman = zeros(Basis.nx,Basis.nx);
Estinfo(1).PRTS = zeros(Basis.nx,Basis.nx);
Estinfo = repmat(Estinfo(1),1,N);

% -------------------------------------------------------
% Initialise parameter estimate
% -------------------------------------------------------

Thetainfo(1).thetaest = zeros(Constants.d,1);
% Thetainfo(1).thetaest = thetatrue;
Thetainfo(1).vartheta = 1*eye(Constants.d);%1000*eye(Constants.d);
Thetainfo(1).muest = 1;
Thetainfo(1).varmu = 0.1;
Thetainfo(1).precisionest = repmat(1/(sigmaW^2),Basis.nx,1);
Thetainfo(1).varprecision = repmat(1,Basis.nx,1);
Thetainfo(1).Meanprecmat = eye(Basis.nx);
% Thetainfo(1).Meanprecmat = inv(FieldMatrices.W);
Thetainfo(1).kernelpar = 0.001;
Thetainfo(1).varkernelpar = 0.1;
Thetainfo(1).rho = 0.5;
Thetainfo(1).varrho = 0.5;

if Estmu == 'n'
    Thetainfo(1).muest = bias;
    Thetainfo(1).varmu = 0;
end
if Estrho == 'n'
    Thetainfo(1).rho = rho;
    Thetainfo(1).varrho = 0;
end
if Esttheta == 'n'
    Thetainfo(1).thetaest = thetatrue;
    Thetainfo(1).vartheta = 0*eye(Constants.d);
end
if Estkernel == 'n'
    Thetainfo(1).kernelpar = 1;
    Thetainfo(1).varkernelpar = 0;
    FieldMatrices.V(:,:,1) = FieldMatrices.PSIx;
end
if Estprec == 'n'
%     my_sum(i) = 0; for i = 1:N my_sum(i) = length(spikeaccept(i).Coords(:,1)); end
%     Educated_var_guess = var(log(my_sum/s(end).^2));
    Thetainfo(1).Meanprecmat = eye(Basis.nx);
end

Thetainfo = repmat(Thetainfo(1),1,numiters);
Constants.theta = b;

% -------------------------------------------------------
% Point process parameters
% -------------------------------------------------------
beta = 1;
Constants.beta = beta;


%-----------------------------------------------
% Inference
%-----------------------------------------------
%Estimate initial condition
EstIntensity = DiggleContSpace(spikeaccept(1).Coords,Constants);
y = Regress(EstIntensity,Basis,Thetainfo(1).muest ,beta);
FieldMatrices.R = sigmaR^2*eye(Basis.nx);
Estinfo(1) = KalmanFilter(y,Estinfo(1),Constants,FieldMatrices,Basis);

% Run VBEM
for m = 2:200
    Estinfo = VBEMSmoother_full_nonlinear(spikeaccept,Constants,Estinfo,FieldMatrices,Basis,Thetainfo(m-1));
    Thetainfo(m) = VBM(Estinfo,Thetainfo(m-1),FieldMatrices,Constants,Basis,spikeaccept,Estmu,Estprec,Estkernel,Esttheta,Estrho);
%     Thetainfo(m).Meanprecmat = diag(Thetainfo(m).precisionest);
    save('LatestResults')
    %Breaking Condition
    if norm(Thetainfo(m).thetaest - Thetainfo(m-1).thetaest) < 0.005 ...
            && (abs(Thetainfo(m).muest - Thetainfo(m-1).muest) < 0.01) ...
            && max(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) < 1.05 ...
            && min(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) > 0.95 ...
            && norm(Thetainfo(m).kernelpar - Thetainfo(m-1).kernelpar) < 0.005
        break
    end
end

save(filename)
toc
%----------------------------------------------------------



function [Estinfo] = VBEMSmoother_full_nonlinear(spikes,Constants,Initinfo,FieldMatrices,Basis,Thetainfo)

% -------------------------------------------------------------
% Function VBEMFilter_full_nonlinear
% inputs:
% outputs:
% -------------------------------------------------------------

PSIxinv = inv(FieldMatrices.PSIx);
SigmaW = inv(Thetainfo.Meanprecmat);
PSIxinvV = PSIxinv*FieldMatrices.V(:,:,1);

rho = Thetainfo.rho;
varrho = Thetainfo.varrho;
Aest = rho*Thetainfo.kernelpar*PSIxinvV;
beta = Constants.beta;

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
AQinvA = (rho^2 + varrho)*(Thetainfo.kernelpar^2 + Thetainfo.varkernelpar)*...
    PSIxinvV'*Qinv*PSIxinvV;
dt = Constants.dt;
y = reshape(spikes,[],size(spikes,3));
Estinfo = Initinfo;

options = foptions;
options(14) = 2000;
% options(9) = 1;
options(2) = 0.1;
options(3) = 0.1;

for i = 2:Constants.N
    xestprior(:,i) = Aest*xestpost(:,i-1) + Constants.dt*theta;
    Sigmaprior(:,:,i) = Aest*Sigmapost(:,:,i-1)*Aest' + SigmaW;
    
    Sigmatilde = inv(inv(Sigmapost(:,:,i-1)) + AQinvA);
    Sigmastar = inv(Qinv - Qinv*Aest*Sigmatilde*Aest'*Qinv);
    mustar = Sigmastar*(Qinv*Aest*Sigmatilde*(inv(Sigmapost(:,:,i-1))*xestpost(:,i-1) - Aest'*Qinv*theta*dt) + Qinv*theta*dt);
    
    %Scan data for NaNs
    spikecoords = spikes(i).Coords;
    
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.tau1(j),Basis.tau2(j))';
        end
    else phieval = zeros(Basis.nx,size(spikecoords,1));
    end
    
    myint = @(xx)  SingleIntForward(xx,Constants,Basis,Thetainfo);
    myint2 = @(xx)  MultiIntForward(xx,Constants,Basis,Thetainfo);
    f = @(xx) -(sum(mu + beta*phieval'*xx') - myint(xx') - ((xx' - mustar)'/Sigmastar)*(xx' - mustar)./2);
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


Sigmabeta(:,:,i) = 9*eye(Basis.nx);
xbeta(:,end) = xestpost(:,end);

for i = Constants.N-1:-1:1
    
    spikecoords = spikes(i+1).Coords;
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.tau1(j),Basis.tau2(j))';
        end
    else phieval = zeros(Basis.nx,size(spikecoords,1));
    end
    
    %     Hacky way linearizing around filtered estimated
    %     Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(xestpost(:,i+1)));
    %     mudash = xestpost(:,i+1) + Sigmadash*(inv(Sigmabeta(:,:,i+1))*(xbeta(:,i+1) - xestpost(:,i+1)) + sum(beta*phieval,2) - myint2(xestpost(:,i+1)));
    
    %   Proper way
    f = @(xx) -(sum(mu + beta*phieval'*xx') - myint(xx') - ((xx' - xbeta(:,i+1))'/Sigmabeta(:,:,i+1))*(xx' - xbeta(:,i+1))./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmabeta(:,:,i+1) + xbeta(:,i+1)'/Sigmabeta(:,:,i+1));
    mudash = scg(f,xestpost(:,i+1)',options,gradf)';
    Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(mudash));
    
    Sigmatilde = inv(inv(Sigmadash) + Qinv);
    Sigmabeta(:,:,i) = inv(AQinvA - Aest'*Qinv*Sigmatilde*Qinv*Aest);
    xbeta(:,i) = Sigmabeta(:,:,i)*(-dt*Aest'*Qinv*theta + Aest'*Qinv*Sigmatilde*(inv(Sigmadash)*mudash + Qinv*dt*theta));
    
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
xestpost(:,1) = Aest*xestRTS(:,2) - theta*dt; %initial state to the next iteration
% Sigmapost(:,:,1) = SigmaW*inv(eye(Basis.nx)-A*A'); %initial state

function [Estinfo] = KalmanFilter(y,Initinfo,Constants,FieldMatrices,Basis)


% -------------------------------------------------------------
% Function KalmanFilter
% inputs:
% outputs:
% Description: In conjunction with Diggle estimator to put into account the
% system dynamics
% -------------------------------------------------------------

SigmaW = 0.2^2*FieldMatrices.W;
SigmaR = FieldMatrices.R;

xestpost = zeros(Basis.nx,Constants.N);
sigma2estpost = zeros(Basis.nx,Basis.nx,Constants.N);
xestprior = zeros(Basis.nx,Constants.N);
sigma2estprior = zeros(Basis.nx,Basis.nx,Constants.N);

sigma2estprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
xestpost(:,1) = Initinfo(1).xestpost;

Estinfo = Initinfo;
C = eye(Basis.nx);

S = C*sigma2estprior(:,:,1)*C' + SigmaR;
K = sigma2estprior(:,:,1)*C'*inv(S);
xestpost(:,1) = xestprior(:,1) + K*(y(:,1) - C*xestprior(:,1));
sigma2estpost(:,:,1) = (eye(Basis.nx) - K*C)*sigma2estprior(:,:,1);
Estinfo(1).xestpost = xestpost(:,1);


function [lambda] = DiggleContSpace(spikes,Constants)

% -------------------------------------------------------------
% Function Diggle
% inputs:
% outputs:
% Description: The Diggle estimator
%---------------------------------------------------------------

lambda = zeros(Constants.J,Constants.J);
frame = lambda;
buff = 1;
rabs = buff*Constants.ds;
% rabs = 3.2;
% buff = rabs/Constants.ds;
r = 4;
for i = 1:size(spikes,1)
    %Find x pixel
    [temp,x] = find((spikes(i,1) - Constants.s).^2 == min((spikes(i,1) - Constants.s).^2));
    %Find y pixel
    [temp,y] = find((spikes(i,2) - Constants.s).^2 == min((spikes(i,2) - Constants.s).^2));
    frame(y,x) = frame(y,x)+1;
end


for i = buff:Constants.J-buff
    for j = buff:Constants.J-buff
        lambda(j,i) = pixincircle(frame,i,j,r)/(pi*rabs^2)/Constants.dt;
    end
end


function y = Regress(lambda,Basis,mu,beta)

phivec = (reshape(Basis.phi,[],Basis.nx))';
y = zeros(Basis.nx,size(lambda,3));
lambda = lambda + 0.001;
for i = 1:size(lambda,3)
    lambdavec = reshape(lambda(:,:,i),[],1);
    y(:,i) = inv(phivec*phivec')*phivec*((log(lambdavec) - mu)/beta);
end

function [Thetainfo] = VBM(Estinfo,Thetainfo,FieldMatrices,Constants,Basis,spikes,Estmu,Estprec,Estkernel,Esttheta,Estrho)

PSIx = FieldMatrices.PSIx;
PSIxinv = inv(PSIx);
SigmaW = FieldMatrices.W;
d = Constants.d;
dt = Constants.dt;
ds = Constants.ds;
N = Constants.N;
PSIxinvV = PSIxinv*FieldMatrices.V(:,:,1);
rho =  Thetainfo(1).rho;
varrho = Thetainfo(1).varrho;
Aest = rho*Thetainfo.kernelpar*PSIxinvV;

%Estimate precision
alpha = 5;
beta = 0.2;

precest = repmat(1/(0.2^2),Basis.nx,1);
varprec = repmat(0,Basis.nx,1);
if Estprec == 'y'
        for j = 1:Basis.nx
            k1 = 0;
            for i = 3:N-2
                k1 = k1+ Estinfo(i).W(j,j) + rho^2*Estinfo(i-1).W(j,j) + dt^2*(Thetainfo(1).vartheta(j,j) + Thetainfo(1).thetaest(j)^2) ...
                    - 2*rho*Estinfo(i-1).S(j,j) + 2*dt*Thetainfo(1).thetaest(j)*Estinfo(i-1).xestRTS(j)*rho ...
                    - 2*dt*Thetainfo(1).thetaest(j)*Estinfo(i).xestRTS(j);
            end
            alpha0 = (Constants.N-4)/2;%alpha0 = alpha + (Constants.N-4)/2;
            beta0 = k1/2;%beta0 = beta + k1/2;
            precest(j) = alpha0/beta0;
            varprec(j) = alpha0/beta0^2;
        end
        Meanprecmat = diag(precest);
    
    %Whole Wishart distribution
    dofprior = 10;
    Prevmatprior = 2/dofprior*eye(Basis.nx);
    Gamma = 0;
    for i = 3:N-2
        Gamma = Gamma + Estinfo(i).W + (rho^2 + varrho)*(Thetainfo(1).kernelpar^2 + Thetainfo(1).varkernelpar)*PSIxinvV*Estinfo(i-1).W*PSIxinvV' ...
            + Thetainfo(1).vartheta + Thetainfo(1).thetaest*Thetainfo(1).thetaest' ...
            - Aest*Estinfo(i-1).S - Estinfo(i-1).S'*Aest' ...
            - Estinfo(i).xestRTS*Thetainfo(1).thetaest' - Thetainfo(1).thetaest* Estinfo(i).xestRTS'  ...
            + Aest*Estinfo(i-1).xestRTS*Thetainfo(1).thetaest' + Thetainfo(1).thetaest*Estinfo(i-1).xestRTS'*Aest';
    end
    Precmatpost = inv(inv(Prevmatprior) + Gamma);
    Meanprecmat = (dofprior + N-4)*Precmatpost;
else
    Meanprecmat = Thetainfo(1).Meanprecmat;
end



%Estimate theta
if Esttheta == 'y'
    v = zeros(d,1);
    for i = 3:Constants.N-2
        v(:,1) = v(:,1) + dt*inv(SigmaW)*(Estinfo(i).xestRTS - Aest*Estinfo(i-1).xestRTS);
    end
    vartheta = SigmaW./dt^2/(Constants.N-4);
    thetaest = vartheta*v;
else
   thetaest = Thetainfo.thetaest;
   vartheta = Thetainfo.vartheta;
end



%Estimaterho
if Estrho == 'y'
    rhoprior = 1;
    varrhoprior = 1;
    v = 0;
    Upsilon = 0;
    for i = 3:N-2
        v = v + trace(Thetainfo(1).Meanprecmat*Estinfo(i-1).S' - Thetainfo(1).Meanprecmat*Thetainfo.thetaest*Estinfo(i-1).xestRTS');
        Upsilon = Upsilon + trace(Thetainfo(1).Meanprecmat*Estinfo(i-1).W);
    end
    Thetainfo(1).varrho = inv(inv(varrhoprior) + Upsilon);
    Thetainfo(1).rho = Thetainfo(1).varrho*(inv(varrhoprior)*rhoprior + v);
end


%EstimateKernel
if Estkernel == 'y'
    kernelprior = 1;
    varkernelprior = 1;
    v = 0;
    Upsilon = 0;
    for i = 3:N-2
        v = v + trace(Thetainfo.Meanprecmat*PSIxinv*FieldMatrices.V(:,:,1)*Estinfo(i).S' - Thetainfo.Meanprecmat*Thetainfo.thetaest*Estinfo(i).xestRTS'*FieldMatrices.V(:,:,1)'*PSIxinv);
        Upsilon = Upsilon + trace(FieldMatrices.V(:,:,1)'*PSIxinv*Thetainfo.Meanprecmat*PSIxinv*FieldMatrices.V(:,:,1)*Estinfo(i).W);
    end
    Thetainfo(1).varkernelpar = inv(inv(varkernelprior) + Upsilon);
    Thetainfo(1).kernelpar = Thetainfo(1).varkernelpar*(inv(varkernelprior)*kernelprior + v);
end

%Estimate mu
if Estmu == 'y'
    muprior = 0;
    Sigmamuprior = 10;
    totalspikenum = 0;
    totalvoidsum = 0;
    phitemp = zeros(Constants.J,Constants.J,Basis.nx);
    Basis_vec_mat = reshape(Basis.phi,[],Basis.nx);
    for i = 3:N-2
        Umean = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS,1,1,Basis.nx)),3);
        MySigma = Estinfo(i).PRTS; %CHANGE TO FULL COVARIANCE
        phitemp_vec = Basis_vec_mat*MySigma;
        phitemp = reshape(phitemp_vec,Constants.J,Constants.J,Basis.nx);
        %         for j = 1:Basis.nx
        %             phitemp(:,:,j) = sum(multiprod(Basis.phi,reshape(MySigma(:,j),1,1,Basis.nx)),3);
        %         end
        phisigma = sum(phitemp.*Basis.phi,3);
        totalvoidsum = totalvoidsum + dt*Constants.ds^2*sum(sum(exp(Constants.beta*Umean + Constants.beta^2*phisigma./2)));
        totalspikenum = totalspikenum + size(spikes(i).Coords,1);
    end
    fmu = @(mu) -(-(mu - muprior)^2/(2*Sigmamuprior) + mu*totalspikenum - exp(mu)*totalvoidsum);
    gradfmu = @(mu) -(-mu/Sigmamuprior + totalspikenum - exp(mu)*totalvoidsum);
    options = foptions;
    options(9) = 1;
    Thetainfo(1).muest = scg(fmu,muprior,options,gradfmu)';
    Thetainfo(1).varmu = 1/(1/Sigmamuprior + exp(Thetainfo(1).muest)*totalvoidsum);
end

Thetainfo(1).precisionest = precest;
Thetainfo(1).varprecision = varprec;
Thetainfo(1).Meanprecmat = Meanprecmat;
Thetainfo.thetaest = thetaest;
Thetainfo.vartheta = vartheta;
%---------------------------
%Auxiliary functions below
%---------------------------
function y = pixincircle(A,xc,yc,radius)

% Engine
y=1:size(A,1);
x=1:size(A,2);
[X Y]=meshgrid(x,y);

% This assume the circle falls *entirely* inside the image
R2 = (X-xc).^2+(Y-yc).^2;
% c = contourc(x,y,R2,[0 0]+radius^2);
[c1,c2] = find(R2 < radius^2);
% c = round(c(:,2:end)); % pixels located ~ on circle
c(1,:) = c1;
c(2,:) = c2;
c = round(c(:,2:end)); % pixels located in circle
Ac = A(sub2ind(size(A),c(1,:),c(2,:))); % extract value
y = sum(Ac);


function [phi] = LocalisedKernelPhi(s1,s2,mu1,mu2,tau1,tau2)
s1 = s1(1,:);
s2 = s2(:,1);
J = size(s1,2);
nx = length(mu1);
phi = zeros(J,J,nx);
for i = 1:nx
    beta1 =tau1(i);
    beta2 = tau2(i);
    l1 = 2*pi/beta1;
    l2 = 2*pi/beta2;
    [temp,ilow1] = min(abs(s1 - (mu1(i) - l1)));  %Find centre of kernel
    [temp,ihigh1] = min(abs(s1 - (mu1(i) + l1)));  %Find centre of kernel
    [temp,ilow2] = min(abs(s2 - (mu2(i) - l2)));  %Find centre of kernel
    [temp,ihigh2] = min(abs(s2 - (mu2(i) + l2)));  %Find centre of kernel
    s_on1 = s1(ilow1:ihigh1);
    s_on2 = s2(ilow2:ihigh2);
    [Delta1,Delta2] = meshgrid(beta1*abs(s_on1 - mu1(i)),beta2*abs(s_on2 - mu2(i)));
    phi(ilow2:ihigh2,ilow1:ihigh1,i) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
end

function [phi] = LocalisedKernelPhi_Cont(s1,s2,mu1,mu2,tau1,tau2)

phi = zeros(length(s1),1);

beta1 = tau1;
beta2 = tau2;
l1 = 2*pi/beta1;
l2 = 2*pi/beta2;
slow1 = (mu1 - l1);  %Find centre of kernel
shigh1 = (mu1 + l1);  %Find centre of kernel
slow2 = (mu2 - l2);  %Find centre of kernel
shigh2 = (mu2 + l2);  %Find centre of kernel

Delta1 = beta1*abs(s1 - mu1);
Delta2 = beta2*abs(s2 - mu2);
phi(:,1) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
phi = phi.*~(s1 < slow1 | s1 > shigh1 | s2 < slow2 | s2 > shigh2);




function y = gaussinner(s1,s2,phi)
nx = size(phi,3);
y = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        y(i,j) = trapz(s2(:,1),(trapz(s1(1,:),phi(:,:,i).*phi(:,:,j),2)));
    end
end

function y = gaussinner2(s1,s2,phi1,phi2)
nx = size(phi1,3);
y = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        y(i,j) = trapz(s2(:,1),(trapz(s1(1,:),phi1(:,:,i).*phi2(:,:,j),2)));
    end
end

function V = FindV(s1,s2,Basis,KernelBasis)
nx = Basis.nx;
d = KernelBasis.d;
V = zeros(d,nx,nx);
for i = 1:d
    for j = 1:nx
        for k = 1:nx
            V(j,k,i) = trapz(s2(:,1),(trapz(s1(1,:),Basis.phi(:,:,j).*KernelBasis.PHI(:,:,i,k),2)));
        end
    end
end




function y = SingleIntForward(x,Constants,Basis,Thetainfo)
beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y = Constants.dt*Constants.ds^2*exp(mu + varmu/2)*sum(exp(beta*Uest));


function y = MultiIntForward2(x,Constants,Basis,Thetainfo)
beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y  = Constants.dt*Constants.ds^2*exp(mu+varmu/2)*beta^2*Basis.Basisvec'*(Basis.Basisvec.*repmat(exp(beta*Uest),1,Basis.nx));


function y = MultiIntForward(x,Constants,Basis,Thetainfo)

beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y =  Constants.dt*Constants.ds^2*exp(mu + varmu/2)*beta*Basis.Basisvec'*exp(beta*Uest);


