function cv = cov_vector_speed(x, V, theta,tau2)
% x - new point (d x 1), d dimension
% V - design points (d x k), d dimension, k points
% theta - factor of Gaussian covariance (d x 1)
% tau2 - variance of radom field

n = size(V,2);
% Gaussian Covariance
cv = repmat(tau2,1,n) .* exp(-sum(repmat(theta,1,n) .* (repmat(x,1,n) - V).^2,1));
   
end