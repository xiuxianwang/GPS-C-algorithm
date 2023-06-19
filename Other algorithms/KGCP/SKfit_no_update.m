function model = SKfit_no_update(X, Y, B, model_old, gammaP)
% fit a stochastic kriging model to simulation output

tau2 =  model_old.tausquared;
beta = model_old.beta;
theta = model_old.theta;
sigma2 = model_old.sigma2;

[k d] = size(X);

if not(all(size(Y)==[k 1]))
    error('Output vector and design point matrix must have the same number of rows.')
end
if (size(B,1)~=k)
    error('Basis function and design point matrices must have the same number of rows.')
end
if not(all(B(:,1)==1))
    error('The first column of the basis function matrix must be ones.')
end
if ((gammaP ~= 1)&(gammaP ~= 2)&(gammaP ~= 3))
    error('please type of correlation function: 1 = exponential, 2 = gauss, 3 = cubic.')
end

% Normalize data by scaling each dimension from 0 to 1
minX = min(X);  
maxX = max(X);
X = (X - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

% calculate the distance matrix between points (copied from DACE)
% distances are recorded for each dimension separately
ndistX = k*(k-1) / 2;        % max number of non-zero distances
ijdistX = zeros(ndistX, 2);  % initialize matrix with indices
distX = zeros(ndistX, d);    % initialize matrix with distances
temp = 0;
for i = 1 : k-1
    temp = temp(end) + (1 : k-i);
    ijdistX(temp,:) = [repmat(i, k-i, 1) (i+1 : k)']; 
    distX(temp,:) = repmat(X(i,:), k-i, 1) - X(i+1:k,:); 
end
IdistX = sub2ind([k k],ijdistX(:,1),ijdistX(:,2));

% distance matrix raised to the power of gamma
D = zeros(k,k,d);
for p=1:d
    temp = zeros(k);
    if (gammaP == 3)
        temp(IdistX) = abs(distX(:,p));
    else
        temp(IdistX) = -abs(distX(:,p)).^gammaP;
    end
    D(:,:,p) = temp+temp';
end


% calculate estimates of the correlation and covariance matrices
if gammaP == 3
    Rhat = corrcubR(theta,D);
else
    Rhat = correxpR(theta,D);
end
%-------------------------------------------------------------------------------%
%estimate tau2, sigma2 using g
Sigmahat  = tau2*Rhat + diag(sigma2*ones(k,1));
Lhat = chol(Sigmahat)';
Lhatinv = inv(Lhat);
Sigmahatinv = Lhatinv'*Lhatinv;

%-------------------------------------------------------------------------------%

% output MLEs and other things useful in prediction
model.tausquared =  tau2;
model.beta = beta;
model.theta = theta;
model.sigma2 = sigma2;
model.X = X;
model.minX = minX;
model.maxX = maxX;
model.gamma = gammaP;
model.Sigma = Sigmahat;
model.L = Lhat;
model.Z = Lhat\(Y-B*beta);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output some additional parameters, added by SHEN Haihui May 23 2017 %%%
% model.negLL = logPL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP);
model.Sigmainv = Sigmahatinv;
model.B = B;
model.Y = Y;
