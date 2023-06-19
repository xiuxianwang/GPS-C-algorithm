function f= SKpredict_fmincon(model,Xpred,A,Bpred,Z,sigma)
% this is objective function that is used in the solver "fmincon" to find
% the maximum solution of the Gaussian-surroagte model

% retrieve model parameters from model structure obtained from SKfit
minX = model.minX;
maxX = model.maxX;
[k d] = size(A);
theta = model.theta;
gammaP = model.gamma;
beta = model.beta;
tau2 = model.tausquared;

% normalization for the design points
A = (A - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

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
        - repmat(A,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(A,[1 1 K]));
end

% calculate correlations between prediction points and design points
D1 = distXpred;
if gammaP == 3
    T = repmat(reshape(theta,[1 d 1]),[k1 1 K]);
    Rpred = tau2*prod(((D1<=(T./2)).*(1-6*(D1./T).^2+6*(D1./T).^3) ...
        +((T./2)<D1 & D1<=T).*(2*(1-D1./T).^3)),2);
else
    Rpred = tau2*exp(sum(-D1.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end
Rpred = reshape(Rpred,[k K 1]);

% calculate responses at prediction points
B1 = ones(size(A,1)^1,1);
f = Bpred*beta + Rpred'*sigma*(Z-B1*beta);
f = -f;

