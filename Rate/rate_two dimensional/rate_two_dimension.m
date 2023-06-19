%% basic parameters of the algorithm
d=2;       %dimension of the problem
MaxSteps=600;      %Maximum iteration of the algorithm
s=10;       %number of solutions sampled in each iteration
lambda = 0.25;  %assign a value for simualtion variance


% determine the random seed
s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

%% generate sample paths of a two-domensional Gaussian process
mu_0 = ones(1,400); % unconditonal mean
sigma2 = 4;    % unconditional variance
V = 0.1*ones(400,1);   % add some noise to avoid numerical issues
% creat a uniform grid with 400 points in the space
A_grid = zeros(400,2);
count = 0;
for i=1:20
    for j = 1:20
        count = count+1;
        A_grid(count,1) = 0.05*i;
        A_grid(count,2) = 0.05*j;
    end
end
% calculate the covariance between the points of the uniform grid
cov = zeros(400,400);
for i =1:size(A_grid,1)
    for j = 1:size(A_grid,1)
        cov(i,j) = sigma2 * exp(-80*((A_grid(i,1)-A_grid(j,1))^2+(A_grid(i,2)-A_grid(j,2))^2));
    end
end
Y = mvnrnd(mu_0,cov);   % generate observarions of the uniform grid 
B = ones(size(A_grid,1)^1,1);
skriging_model = SKfit_create_sample_path(A_grid,Y',B, V, 2);   %fit a Gaussian-surrogate model based on the observarions

% identify the maximum value of the generated sample path
% construct a dense grid for evaluation
A_eva = zeros(101*101,2);
count = 0;
for i=1:101
    for j = 1:101
        count = count+1;
        A_eva(count,1) = 0.01*(i-1);
        A_eva(count,2) = 0.01*(j-1);
    end
end
B_eva = ones(size(A_eva,1)^1,1);
Y_true = SKpredict_sample_path(skriging_model,A_eva,B_eva);  % evaluate the dense grid
[true_max_value,location] = max(Y_true);  % use the maximum value of the dense grid as the global optimal value of the generates sample path
 
%% The GPS-C algorithm
% generate vectors to store     
Opti_Record = zeros(1,MaxSteps);
Opti_Record_gap = zeros(1,MaxSteps);

% Step 0 - initializaiton, assign values for the Gaussian process parameters
Sk = [0,0];      % aritrarily sample a point
Gk = [0];
% fit the initial kriging model
B = ones(size(Sk,1)^1,1);       
skriging_model_2 = SKfit(Sk,Gk,1,lambda,2);   %fit a kriging model


% Step 1 - sample some design points following a uniform distribution 
Sample = rand(s,2);

% take observations
SSample = [];
G = [];
G = fun_many_loc(Sample,skriging_model);
SSample = Sample;
SSample(:,d+1) = ones(size(Sample,1),1);
SSample(:,d+2) = G;
SSample(:,d+3) = ones(size(Sample,1),1)*lambda;

% matrix inversion
sigma = Matrix_inverse(skriging_model_2,SSample(:,d+3),SSample(:,1:d));

% use the dense grid to find the current best solution of each iteration
A = A_eva;
B = ones(size(A,1)^1,1);
[A(:,d+1),sigma]  = SKpredict(skriging_model_2,A(:,1:d),Sample,B,SSample(:,d+2),sigma);
[g_star,max_location] = max(A(:,d+1));
x_star = A(max_location,1:d);

step = 1;
Opti_Record(1) = g_star;
Opti_Record_gap(1) = abs(g_star-true_max_value);
    
% Step2 - iteration
while step<MaxSteps

    % sample s points from the constructed distribution
    Sample = Sampling(skriging_model_2,SSample,g_star,x_star,s,sigma);
    
    % take one observation for each newly sampled point
    G = [];
    G = fun_many_loc(Sample,skriging_model);

    % Merge to SSample
    for i = 1:size(Sample,1)
        m = size(SSample,1);
        SSample(m+1,d+2) = G(i);
        SSample(m+1,1:d+1) = [Sample(i,:),1]; 
        SSample(m+1,d+3) = lambda;
    end

    % matrix inversion
    sigma = Matrix_inverse(skriging_model_2,SSample(:,d+3),SSample(:,1:d));

    % use the dense grid to find the current best solution of each iteration
    B = ones(size(A,1)^1,1);
    [A(:,d+1),sigma]  = SKpredict(skriging_model_2,A(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
    [g_star,max_location] = max(A(:,d+1));
    x_star = A(max_location,1:d);

    step = step+1;
    fprintf('%d  [%.2f, %.2f] %.2f \n',step,x_star,g_star-true_max_value); % output current best solution

    % record the result of each step
    Opti_Record(step) = g_star;
    Opti_Record_gap(step) = abs(g_star-true_max_value);
end


