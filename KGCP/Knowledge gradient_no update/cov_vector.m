function cv = cov_vector(x_sample, sample_x, model)

K = size(x_sample,1);     % number of prediction points
% calculate distance matrix for prediction points
minX = model.minX;
maxX = model.maxX;
x_sample = (x_sample - repmat(minX,K,1)) ./ repmat(maxX-minX,K,1);
K1 = size(sample_x,1);     % number of prediction points
sample_x = (sample_x - repmat(minX,K1,1)) ./ repmat(maxX-minX,K1,1);

x = x_sample';
V = sample_x';
theta = model.theta;
tau2 = model.tausquared;
n = size(V,2);
% Gaussian Covariance
cv = repmat(tau2,1,n) .* exp(-sum(repmat(theta,1,n) .* (repmat(x,1,n) - V).^2,1));
   
end