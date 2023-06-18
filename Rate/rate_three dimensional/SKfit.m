function model = SKfit(X, Y, B, Vhat, gammaP)
% all the parameters are given 

% output MLEs and other things useful in prediction
model.tausquared =  9;
model.beta = 1;
model.theta = [40;40;40];
model.X = X;
model.minX = [0,0,0];
model.maxX = [1,1,1];
model.gamma = 2;