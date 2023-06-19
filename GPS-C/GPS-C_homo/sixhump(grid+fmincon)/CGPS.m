%%%% Continuous Gaussian Process-Based Search(GPS-C) Algorithm
% This is the code of the GPS-C algorithm, which is an adaptive random search based simulation optimization
% algorithm for continuous problems. Since the Gaussian-surrogate model is
% used in this algorithm for 1) estimating the unknown objective function,
% and 2) constructing sampling distribution, the code of Stochastic kriging
% of Ankenman et. al. (2010) is used.
clear,clc, close all;

% setting of the problem
fun_name = 'Fun_Sixhump';
 % get the information of the test problem
switch fun_name
    case 'Fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=1.031628;
    case 'Fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'Fun_Hills'
        num_vari=2; design_space=[0,0;100,100];                  optimum=20;
    case 'Fun_Rosenbrock'
        num_vari=10; design_space=[-10,-10,-10,-10,-10,-10,-10,-10,-10,-10;10,10,10,10,10,10,10,10,10,10];            optimum=0;
    otherwise
        error('objective function is not defined!')
end  
%% basic parameters of the algorithm

d=2;       %dimension of the problem
MaxSteps=78;      %Maximum iteration of the algorithm
s=10;       %number of solutions sampled in each iteration
% generate uniform grids in the feasible region for evaluation
count = 0;
A = zeros(51*51,2);
for i=1:51   % the grid for the sixhump problem
    for j = 1:51
        count = count+1;
        A(count,1) = -2+0.08*(i-1);
        A(count,2) = -2+0.08*(j-1);
    end
end

%% The main steps of the GPS-C algorithm

s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

% Step 0 - initializaiton
Opti_Record = zeros(MaxSteps,1); % record the current best solution of each iteration
SSample = []; % record the information of each sampled design point
G = []; 

% initial design points using Latin hypercube sampling method
num_initial = 20;
Sample_ini = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
G_ini = feval(fun_name,Sample_ini);

% randomly sample s points in the feasible region (optional)
Sample = 4*(rand(s,d)-0.5*ones(s,d));
G = feval(fun_name,Sample); % take observations of the randomly sampled points

% record the observations 
SSample = [Sample_ini;Sample];
SSample(:,d+1) = ones(size(SSample,1),1); % the number of observations of each point
SSample(:,d+2) = [G_ini;G];

% fit a krging model
B = ones(size(SSample,1)^1,1);       % basis function matrix at design points
skriging_model = SKfit(SSample(:,1:d),SSample(:,d+2),B,2);   %fit a kriging model
lambda = skriging_model.sigma2;  % the estimated variance of simulation noise
SSample(:,d+3) = ones(size(SSample,1),1)*lambda;

% matrix inverse for further calculating the conditional mean and variance
sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d)); 

% find current best solution using grid search + fmincon
B = ones(size(A,1)^1,1);
Y_hat = SKpredict(skriging_model,A(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
[start_value,max_location] = max(Y_hat); % find the optimal solution of the grid
start = A(max_location,1:d);  % use the grid-best solution as the start point
myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
lb = [max(start(1,1) - 0.08,-2),max(start(1,2)-0.08,2)];
ub = [min(start(1,1) + 0.08,2),min(start(1,2)+0.08,2)];
B = ones(1^1,1);
x_star = fmincon(@(x) SKpredict_fmincon(skriging_model,x,SSample(:,1:d),B,SSample(:,d+2),sigma),...
                     start,[],[],[],[],lb,ub,[],myopt);
g_star = Fun_Sixhump_free(x_star); % calulate the true objective function value

% Step2 - iteration
step = 1;
Opti_Record(1) = g_star;

while step<MaxSteps
    % sample s points from fitted distribution
    Sample = Sampling(skriging_model,SSample,g_star,x_star,s,sigma);
    G = [];
    G = feval(fun_name,Sample);
    
    % Merge to SSample
    for i = 1:size(Sample,1)
            m = size(SSample,1);
            SSample(m+1,d+2) = G(i);
            SSample(m+1,1:d+1) = [Sample(i,:),1];
    end
    % continue updating the Gaussian process parameters in early stage
    if (step<=10)
        B = ones(size(SSample,1)^1,1);      
        skriging_model = SKfit(SSample(:,1:d),SSample(:,d+2),B,2);  
        lambda = skriging_model.sigma2;
    end
    SSample(:,d+3) = lambda*ones(size(SSample,1),1);
    
    % matrix inversion
    sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d));
    
    % find current best solution using grid search + fmincon
    B = ones(size(A,1)^1,1);
    Y_hat = SKpredict(skriging_model,A(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
    [start_value,max_location] = max(Y_hat);
    start = A(max_location,1:d);
    myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
    lb = [max(start(1,1) - 0.08,-2),max(start(1,2)-0.08,-2)];
    ub = [min(start(1,1) + 0.08,2),min(start(1,2)+0.08,2)];
    B = ones(1^1,1);
    x_star = fmincon(@(x) SKpredict_fmincon(skriging_model,x,SSample(:,1:d),B,SSample(:,d+2),sigma),...
                         start,[],[],[],[],lb,ub,[],myopt);                 
    g_star = Fun_Sixhump_free(x_star);
    
    step = step+1; % update the iteration count
    Opti_Record(step) = g_star;  % record the objective value of current best solution
    fprintf('%d  [%.3f, %.3f] %f \n',step,x_star(1),x_star(2),g_star); % output the current best solution

end
