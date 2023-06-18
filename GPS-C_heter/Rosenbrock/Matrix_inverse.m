function sigma = Matrix_inverse(model,Vhat,A)

% retrieve model parameters from model structure obtained from SKfit
minX = model.minX;
maxX = model.maxX;
[k d] = size(A);
theta = model.theta;
gammaP = model.gamma;
beta = model.beta;
tau2 = model.tausquared;

% normalization of the design points
A = (A - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

% calculate the distance matrix between points yuege(copied from DACE)
% distances are recorded for each dimension separately
ndistA = k*(k-1) / 2;        % max number of non-zero distances
ijdistA = zeros(ndistA, 2);  % initialize matrix with indices
distA = zeros(ndistA, d);    % initialize matrix with distances
temp = 0;
for i = 1 : k-1
    temp = temp(end) + (1 : k-i);
    ijdistA(temp,:) = [repmat(i, k-i, 1) (i+1 : k)']; 
    distA(temp,:) = repmat(A(i,:), k-i, 1) - A(i+1:k,:); 
end
IdistA = sub2ind([k k],ijdistA(:,1),ijdistA(:,2));

% distance matrix raised to the power of gamma
D = zeros(k,k,d);
for p=1:d
    temp = zeros(k);
    if (gammaP == 3)
        temp(IdistA) = abs(distA(:,p));
    else
        temp(IdistA) = -abs(distA(:,p)).^gammaP;
    end
    D(:,:,p) = temp+temp';
end

V = diag(Vhat);
Rhat = correxpR(theta,D);
Sigmahat = tau2*Rhat + V;

% calculate inverse matrix
sigma = inv(Sigmahat);

end