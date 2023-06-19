function model = SKfit2(X, Y, B, Vhat, gammaP)
% the same as the SKfit function except that the estimated conditional
% variance is increased to explore the feasible region better in the early 
% stage of the GPC-S algorithm.
% Note that this funciton is only used in the early stage of the algorithm

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
if (not(all(size(Vhat)==[k 1])))
    error('The variance vector and design point matrix must have the same number of rows.')
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

% diagonal intrinsic variance matrix
V = diag(Vhat);

% initialize parameters theta, tau2 and beta
% inital extrinsic variance = variance of ordinary regression residuals
betahat = (B'*B)\(B'*Y);
tau2_0 = var(Y-B*betahat);

% see stochastic kriging tutorial for explanation
% make correlation = 1/2 at average distance
average_distance = mean(abs(distX));
theta_0 = zeros(d,1);
if gammaP == 3
    theta_0(1:d) = 0.5;
else
    theta_0(1:d) = (log(2)/d)*(average_distance.^(-gammaP));
end

% lower bounds for parameters tau2, theta, and beta included 
% only to satisfy fmincon  
lbtau2 = 0.001*tau2_0;     % naturally 0 
lbtheta = 100*ones(d,1); % naturally 0; increase to avoid numerical trouble
lb = [lbtau2;lbtheta];
ub = [100;300*ones(d,1)];

% maximize profile log-likelihood function (-"logPL")
% subject to lower bounds on tau2 and theta
% myopt = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',500);
myopt = optimset('Display','off','MaxFunEvals',1000000,'MaxIter',500);
parms = fmincon(@(x) logPL(x,k,d,D,B,V,Y,gammaP),...
        [tau2_0;theta_0],[],[],[],[],lb,ub,[],myopt); 

% record MLEs for tau2 and theta 
tau2hat = parms(1);
thetahat = parms(2:length(parms));

% MLE of beta is known in closed form given values for tau2, theta
% and is computed below
% calculate estimates of the correlation and covariance matrices
if gammaP == 3
    Rhat = corrcubR(thetahat,D);
else
    Rhat = correxpR(thetahat,D);
end
Sigmahat  = tau2hat*Rhat + V;
Lhat = chol(Sigmahat)';
Lhatinv = inv(Lhat);
Sigmahatinv = Lhatinv'*Lhatinv;
betahat = inv(B'*Sigmahatinv*B)*(B'*(Sigmahatinv*Y));

% issue warnings related to constraints 
warningtolerance = 0.001;
if min(abs(lbtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial lower bound');
end
if abs(lbtau2 - tau2hat) < warningtolerance
    warning('tau2hat was very close to artificial lower bound');
end

% output MLEs and other things useful in prediction
model.tausquared =  tau2hat*4;  % the estimated tau2 is enlarged for better exploration in early stage of the algorithm (In Hills problem, the tau2 may be underestimated with small amount of design points)
model.beta = betahat;
model.theta = thetahat;
model.X = X;
model.minX = minX;
model.maxX = maxX;
model.gamma = gammaP;
model.L = Lhat;
model.Z = Lhat\(Y-B*betahat);
