function [f, MSE] = SKpredict(model,Xpred,Bpred)
% make predictions at prediction points using a stochastic kriging model  

% retrieve model parameters from model structure obtained from SKfit
X = model.X;
minX = model.minX;
maxX = model.maxX;
[k d] = size(X);
theta = model.theta;
gammaP = model.gamma;
sigma2 = model.sigma2;
beta = model.beta;
Z = model.Z;
L = model.L;
tau2 = model.tausquared;
Sigma = model.Sigma;
Sigmainv = model.Sigmainv;
Y = model.Y;
B = model.B;

% simple check for dimensions of Xpred and X
K = size(Xpred,1);     % number of prediction points
if (size(Xpred,2)~=d)
    error('Prediction points and design points must have the same dimension (number of columns).');
end
if (size(Bpred,1)~=K)
    error('Basis function and prediction point matrices must have the same number of rows.')
end
if not(all(Bpred(:,1)==1))
    error('The first column of the basis function matrix must be ones.')
end

% calculate distance matrix for prediction points
Xpred = (Xpred - repmat(minX,K,1)) ./ repmat(maxX-minX,K,1);
if gammaP == 2
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K]));
end

% calculate correlations between prediction points and design points
D = distXpred;
if gammaP == 3
    T = repmat(reshape(theta,[1 d 1]),[k 1 K]);
    Rpred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
else
    Rpred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
Rpred = reshape(Rpred,[k K 1]);

% calculate responses at prediction points 
f = Bpred*beta + Rpred'*(L'\Z);

if K > 1
    kk =size(Sigma,1);
    BB = ones(kk,1);
    M_temp = [0,BB';BB,Sigma];
    MSE=zeros(K,1);
    for i=1:1:K
        MSE(i,1)=tau2-[1,Rpred(:,i)']*(M_temp\[1;Rpred(:,i)]);
    end
else
    kk =size(Sigma,1);
    BB = ones(kk,1);
    M_temp = [0,BB';BB,Sigma];
    MSE = tau2-[1,Rpred']*(M_temp\[1;Rpred]);
end
