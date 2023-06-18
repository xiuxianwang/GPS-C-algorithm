function f = logPL(parms,k,d,Dist,B,Y,gammaP)

% negative of log profile likelihood  
% This is a function of theta and tau2.
% we have profiled out over beta
% parms = [tausq;theta] 
% k = number of design points
% d = dimension of space
% Dist = matrix of distances between design points
% B = basis functions
% V = intrinsic covariance matrix
% Y = simulation output

if ((length(parms)~=d+1))
    error('Parameter vector must be length d+1');
end
if not(all(size(Y)==[k 1]))
    error('Output vector must be k by 1.')
end
if (size(B,1)~=k)
    error('Basis function matrix must have k rows.')
end
if not(all(B(:,1)==1))
    error('The first column of the basis function matrix must be ones.')
end

if(parms(1)<=0 || min(parms(2:d+1)) <= 0.001)
    f = inf;
    return;
end
%-------------------------------------------------------------------------------------%
% % estimate tau2, sigma2, theta directly
% tausq = parms(1);
% sigmasq = parms(2);
% theta = parms(3:d+2);
% 
% % get correlation matrix given theta
% if gammaP == 3
%     R = corrcubR(theta,Dist);
% else
%     R = correxpR(theta,Dist);
% end
%     
% % sum of extrinsic and intrinsic covariances
% V  = tausq.*R + diag(sigmasq*ones(k,1));
% Rhat = V./(tausq+sigmasq);
% [U,pd] = chol(V);
% [U2,pd2] = chol(Rhat);
% 
% if(pd>0||pd2>0)
%     save data;
%     error('covariance matrix is nearly singular');
% end
% 
% % invert V via Cholesky factorization
% L1 = U';
% L1inv = inv(L1);
% SinV = L1inv'*L1inv;
% 
% %invert Rhat via Cholesky factorization
% L2 = U2';
% L2inv = inv(L2);
% SinRhat = L2inv'*L2inv;
% 
% % the optimal beta given theta and tau2
% beta = inv(B'*SinV*B)*(B'*(SinV*Y)); 
% Z = L2\(Y-B*beta);
% 
% %negative log likelihood
% f =  1/k * log(det(Rhat)) + log(1/k * ((Y-B*beta)'*SinRhat*(Y-B*beta)));
%------------------------------------------------------------------------------------------------------%
% negative log likelihood
%f = (log(det(L)) + 0.5*Z'*Z + 0.5*k*log(2*pi));

%------------------------------------------------------------------------------------------------------%
%estimate g first
g = parms(1);
theta = parms(2:d+1);
% get correlation matrix given theta
if gammaP == 3
    R = corrcubR(theta,Dist);
else
    R = correxpR(theta,Dist);
end
Rhat = R.*g;
for i = 1:size(R,1)
    Rhat(i,i) = 1;
end
beta = inv(B'*(Rhat\B))*(B'*(Rhat\Y));

f = 1/k * log(det(Rhat)) + log(1/k * ((Y-B*beta)'*(Rhat\(Y-B*beta))));
%------------------------------------------------------------------------------------------------------%

