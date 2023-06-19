function sigma = Matrix_inverse_block(model,Vhat,A,sigma_old)

% retrieve model parameters from model structure obtained from SKfit
minX = model.minX;
maxX = model.maxX;
[k, d] = size(A);
theta = model.theta;
gammaP = model.gamma;
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

% block matrix inversion
k1 = size(sigma_old,1);
M_B = Sigmahat(1:k1,k1+1:k);
M_C = Sigmahat(k1+1:k,1:k1);
M_D = Sigmahat(k1+1:k,k1+1:k);
S_D = inv(M_D - M_C * sigma_old * M_B);
temp_M_1 = M_C * sigma_old;
S_B = - sigma_old * M_B * S_D;
S_A = sigma_old + (-S_B) * temp_M_1;
S_C = - S_D * temp_M_1;
sigma = [S_A S_B;S_C S_D];

end