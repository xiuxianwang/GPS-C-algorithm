% calculate the conditional variance 
function y = CalculateMSE(model,x,sigma,A,Bpred)

% retrieve model parameters from model structure obtained from SKfit
X = model.X;
minX = model.minX;
maxX = model.maxX;
[k, d] = size(A);
theta = model.theta;
tau2 = model.tausquared;

% normalization of design points
A = (A - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

% dimensions of prediction points
K=size(x,1);

% calculate distance matrix for prediction points
xpred=(x-repmat(minX,K,1))./repmat(maxX-minX,K,1);
distxpred=abs(repmat(reshape(xpred',[1 d K]),[k,1,1])-repmat(A,[1 1 K])).^2;

% calculate correlations between prediction points and design points
Rpred=tau2*exp(sum(-distxpred.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
Rpred=reshape(Rpred,[k K 1]);

% calculate the conditional variance
if K > 1
    y=zeros(K,1);
    for i=1:1:K
        y(i)=tau2-Rpred(:,i)'*sigma*Rpred(:,i);
    end
    y=y';    
else
    y = tau2-Rpred'*sigma*Rpred;
end

end

