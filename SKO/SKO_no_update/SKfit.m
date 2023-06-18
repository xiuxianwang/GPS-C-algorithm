function model = SKfit(X, Y, B, gammaP)
% fit a stochastic kriging model to simulation output

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

% parameters for sixhump problem
% theta_0 = [3;4];
% g_0 = 0.5;
% lbtheta = [1;1]; 
% lb = [0.01;lbtheta];
% ub = [1;[10;10]];
% parameters for Hills problem
theta_0 = [200;200];
g_0 = 0.5;
lbtheta = [100;100]; 
lb = [0.01;lbtheta];
ub = [1;[300;300]];
% maximize profile log-likelihood function (-"logPL")
% subject to lower bounds on tau2 and theta
myopt = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',500);
parms = fmincon(@(x) logPL(x,k,d,D,B,Y,gammaP),...
        [g_0;theta_0],[],[],[],[],lb,ub,[],myopt); 

g = parms(1);
thetahat = parms(2:length(parms));


% MLE of beta is known in closed form given values for tau2, theta
% and is computed below
% calculate estimates of the correlation and covariance matrices
if gammaP == 3
    Rhat = corrcubR(thetahat,D);
else
    Rhat = correxpR(thetahat,D);
end
%-------------------------------------------------------------------------------%
%estimate tau2, sigma2 using g
RR = g.*Rhat;
for i = 1:size(RR,1)
    RR(i,i) = 1;
end
betahat = inv(B'*(RR\B))*(B'*(RR\Y));
sigmahat2 = 1/k * (Y - B*betahat)'*(RR\ (Y - B*betahat));
tau2hat = sigmahat2 * g;
sigma2hat = sigmahat2 * (1-g);
Sigmahat  = tau2hat*Rhat + diag(sigma2hat*ones(k,1));
Lhat = chol(Sigmahat)';
Lhatinv = inv(Lhat);
Sigmahatinv = Lhatinv'*Lhatinv;

%-------------------------------------------------------------------------------%


% issue warnings related to constraints 
warningtolerance = 0.001;
if min(abs(lbtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial lower bound');
end

% output MLEs and other things useful in prediction
model.tausquared =  tau2hat;
model.beta = betahat;
model.theta = thetahat;
model.sigma2 = sigma2hat;
model.X = X;
model.minX = minX;
model.maxX = maxX;
model.gamma = gammaP;
model.Sigma = Sigmahat;
model.L = Lhat;
model.Z = Lhat\(Y-B*betahat);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output some additional parameters, added by SHEN Haihui May 23 2017 %%%
% model.negLL = logPL([tau2hat; thetahat],k,d,D,B,V,Y,gammaP);
model.Sigmainv = Sigmahatinv;
model.B = B;
model.Y = Y;