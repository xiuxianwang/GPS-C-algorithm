%%%% Continuous Gaussian Process-Based Search(GPS-C) Algorithm
clear,clc, close all;

% setting of the problem
fun_name = 'Fun_Rosenbrock';
 % get the information of the test problem
switch fun_name
    case 'Fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=-1.031628;
    case 'Fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'Fun_Hills'
        num_vari=2; design_space=[0,0;100,100];                  optimum=20;
    case 'Fun_Rosenbrock'
        num_vari=10; design_space=[-10,-10,-10,-10,-10,-10,-10,-10,-10,-10;10,10,10,10,10,10,10,10,10,10];            optimum=0;
    otherwise
        error('objective function is not defined!')
end  
%--------------------------------------------------------------------------

d=10;       %dimension of the problem
MaxSteps = 320;      %Maximum iteration of the algorithm
r0 = 5;    % initial bandwidth of the kernel
beta = 0.009;    % contracting rate of the kernel bandwidth
s=10;       %number of solutions sampled in each iteration

s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

% Step 0 - initializaiton
Opti_Record = zeros(MaxSteps,1);
Opti_Record_true = zeros(MaxSteps,1);
Opti_Record_Point = zeros(MaxSteps,d);

% initial design points using Latin hypercube sampling method
num_initial = 400;
Sample_ini = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
Sample_ini1 = Sample_ini + 2*(rand(num_initial,d)-0.5*ones(num_initial,d));
Sample_ini = [Sample_ini;Sample_ini1];
G_ini = feval(fun_name,Sample_ini);

% randomly sample s points in the feasible region (optional)
Sample = 20*(rand(s,d)-0.5*ones(s,d));

% take and record one observation of each sampled point
SSample = [];    % record the information of each sampled design point
dist_matrix = [];
G = [];
G = feval(fun_name,Sample);    % take observations of the randomly sampled points
SSample = [Sample_ini;Sample];
SSample(:,d+1) = ones(size(SSample,1),1);    % the number of observations of each point
SSample(:,d+2) = [G_ini;G];

% calculate sample variance using kernel method
k = size(SSample,1);
for i = 1:k
    for j = 1:k
        dist_matrix(i,j) = max(abs(SSample(i,1:d)-SSample(j,1:d)));    % calculate the distance matrix of design points
    end 
end
for i = 1:k
    ind = [];
    ind = find(dist_matrix(i,:)<r0);    % identify the points within the kernel
    ln = size(ind,2);
    m_hat = sum(SSample(ind,d+2))/ln;
    lambda_hat = sum( (SSample(ind,d+2) - m_hat*ones(ln,1)).^2 )/ln;
    if (lambda_hat<0.1)    % set a lower bound for the variance of simulation noise
        SSample(i,d+3) = 0.1;
    elseif (lambda_hat>5)
        SSample(i,d+3) = 5;
    else
         SSample(i,d+3) = lambda_hat;
    end
end

% fit a new model
B = ones(size(SSample,1)^1,1);       % basis function matrix at design points
skriging_model = SKfit(SSample(:,1:d),SSample(:,d+2),B, SSample(:,d+3), 2);   %fit a kriging model

% matrix inversion
sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d));
sigma_old = sigma;

% find current sample best
B = ones(size(SSample,1)^1,1);
SSample(:,d+4) = SKpredict(skriging_model,SSample(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
[g_star_1,max_location] = max(SSample(:,d+4));
x_star_1 = SSample(max_location,1:d);

 
% Step2 - iteration
step = 1;

Opti_Record(1) = g_star_1;
Opti_Record_true(1) = Fun_Rosenbrock_free(x_star_1);    % calulate the true objective function value
Opti_Record_Point(1,1:d) = x_star_1;

while step<MaxSteps
    
    % sample s points from fitted distribution
    tic;
    Sample = Sampling(skriging_model,SSample,g_star_1,x_star_1,s,sigma);
    toc;
    % take r=1 observation for each newly sampled point
    G = [];
    G = feval(fun_name,Sample);
    rk =r0 * (k + 1)^(-beta);    % update the bandwidth of the kernel
    
    % Merge to SSample
    for i = 1:size(Sample,1)
            m = size(SSample,1);
            SSample(m+1,d+2) = G(i);
            SSample(m+1,1:d+1) = [Sample(i,:),1]; 
    end
    k = size(SSample,1);
    tic;
    % update the distance matrix
    for j = 1:k
        for i = 1:s
            dist_matrix(j,k-s+i) = max(abs(SSample(j,1:d)-SSample(k-s+i,1:d)));
            dist_matrix(k-s+i,j) = max(abs(SSample(k-s+i,1:d)-SSample(j,1:d)));
        end
    end
    % update the estimated variance at each design point
    for i = 1:k
        ind = [];
        ind = find(dist_matrix(i,:)<rk);
        ln = size(ind,2);
        m_hat = sum(SSample(ind,d+2))/ln;
        lambda_hat = sum( (SSample(ind,d+2) - m_hat*ones(ln,1)).^2 )/ln;
        if (lambda_hat<0.1)
            SSample(i,d+3) = 0.1;
        elseif (lambda_hat>5)
            SSample(i,d+3) = 5;
        else
            SSample(i,d+3) = lambda_hat;
        end
    end
    toc;
    
    % matrix inversion using block matrix inversion
    sigma = Matrix_inverse_1(skriging_model,SSample(:,d+3),SSample(:,1:d),sigma_old);
    sigma_old = sigma;
    
    % find current current sample best solution
    B = ones(size(SSample,1)^1,1);
    SSample(:,d+4) = SKpredict(skriging_model,SSample(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
    [g_star_sample_best,max_location] = max(SSample(:,d+4));
    x_sample_best = SSample(max_location,1:d);
    g_star_1 = g_star_sample_best;
    x_star_1 = x_sample_best;
    
    step = step+1;
    Opti_Record(step) = g_star_1;
    Opti_Record_true(step) = Fun_Rosenbrock_free(x_star_1);
    Opti_Record_Point(step,1:d) = x_star_1;   
    fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f] %.4f \n',step,x_star_1(1),x_star_1(2),x_star_1(3),x_star_1(4),x_star_1(5),x_star_1(6),x_star_1(7),x_star_1(8),x_star_1(9),x_star_1(10),Opti_Record_true(step));

end

