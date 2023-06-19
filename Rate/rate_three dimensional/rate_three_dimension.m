%% basic parameters of the algorithm
d=3;       %dimension of the problem
MaxSteps=600;      %Maximum iteration of the algorithm
s=10;       %number of solutions sampled in each iteration
lambda = 0.25;  %assign a value for simualtion variance


% determine the random seed
s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

%% generate sample paths of a three-domensional Gaussian process
mu_0 = ones(1,1000);   % unconditonal mean
sigma2 = 9;   % unconditional variance
V = 0.2*ones(1000,1);   % add some noise to avoid numerical issues
% creat a uniform grid with 1000 points in the space
A_grid = zeros(1000,3);
for i=1:10
    A_grid(100*(i-1)+1:100*i,1) = 0.1*i;
    for j = 1:10
        A_grid(100*(i-1)+10*(j-1)+1:100*(i-1)+10*j,2) = 0.1*j;
        for k = 1:10
            A_grid(100*(i-1)+10*(j-1)+k,3) = 0.1*k;
        end
    end
end
% calculate the covariance between the points of the uniform grid
cov = zeros(1000,1000);
for i =1:size(A_grid,1)
    for j = 1:size(A_grid,1)
        cov(i,j) = sigma2 * exp(-40*((A_grid(i,1)-A_grid(j,1))^2+(A_grid(i,2)-A_grid(j,2))^2+(A_grid(i,3)-A_grid(j,3))^2));
    end
end
Y = mvnrnd(mu_0,cov);  % generate observarions of the uniform grid 
B = ones(size(A_grid,1)^1,1);
skriging_model = SKfit_create_sample_path(A_grid,Y',B, V, 2);   %fit a Gaussian-surrogate model based on the observarions

% identify the maximum value of the generated sample path
% construct a dense grid for evaluation
A_eva = zeros(101*101*101,3);
for i=1:101
    A_eva(10201*(i-1)+1:10201*i,1) = 0.01*(i-1);
    for j = 1:101
        A_eva(10201*(i-1)+101*(j-1)+1:10201*(i-1)+101*j,2) = 0.01*(j-1);
        for k = 1:101
            A_eva(10201*(i-1)+101*(j-1)+k,3) = 0.01*(k-1);
        end
    end
end
true_max_value = 0;
for i=1:101   % evaluate the dense grid
    B_eva = ones(size(A_eva(10201*(i-1)+1:10201*i,:),1)^1,1);
    Y_true = SKpredict_sample_path(skriging_model,A_eva(10201*(i-1)+1:10201*i,:),B_eva);
    [true_max_value_temp,location] = max(Y_true);
    if true_max_value_temp>true_max_value
        true_max_value = true_max_value_temp;   % use the maximum value of the dense grid as the global optimal value of the generates sample path
    end
end
    
%% The GPS-C algorithm
% generate vectors to store     
Opti_Record = zeros(1,MaxSteps);
Opti_Record_gap = zeros(1,MaxSteps);
    
% Step 0 - initializaiton, assign values for the Gaussian process parameters
Sk = rand(1,3);
Gk = [1];
skriging_model_2 = SKfit(Sk,Gk,1,lambda,2);   %fit a kriging model


% Step 1 - sample some design points following a uniform distribution 
Sample = rand(s,3);
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
sigma_old = sigma;
    
% use the dense grid to find the current best solution of each iteration
g_star = 0;
for i=1:101   % evaluate the dense grid
    AA = A_eva(10201*(i-1)+1:10201*i,:);
    B_eva = ones(size(AA,1)^1,1);
    Y_hat = SKpredict(skriging_model_2,AA,SSample(:,1:d),B_eva,SSample(:,d+2),sigma);
    [g_star_temp,location] = max(Y_hat);
    if g_star_temp>g_star
        g_star = g_star_temp;   % use the maximum value of the dense grid as the global optimal value of the generates sample path
        x_star = AA(location,:);
    end
end

step = 1;
Opti_Record(1) = g_star;
Opti_Record_gap(1) = abs(g_star-true_max_value);

% Step2 - iteration
while step<MaxSteps

    % sample s points from fitted distribution
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
    sigma = Matrix_inverse_block(skriging_model_2,SSample(:,d+3),SSample(:,1:d),sigma_old);
    sigma_old = sigma;

    % use the dense grid to find the current best solution of each iteration
    g_star = 0;
    for i=1:101   % evaluate the dense grid
        AA = A_eva(10201*(i-1)+1:10201*i,:);
        B_eva = ones(size(AA,1)^1,1);
        Y_hat = SKpredict(skriging_model_2,AA,SSample(:,1:d),B_eva,SSample(:,d+2),sigma);
        [g_star_temp,location] = max(Y_hat);
        if g_star_temp>g_star
            g_star = g_star_temp;   % use the maximum value of the dense grid as the global optimal value of the generates sample path
            x_star = AA(location,:);
        end
    end

    step = step+1;
    fprintf('%d  [%.2f, %.2f, %.2f] %.2f \n',step,x_star,g_star-true_max_value);  % output current best solution
    
    % record the result of each step
    Opti_Record(step) = g_star;
    Opti_Record_gap(step) = abs(g_star-true_max_value);
end
    

