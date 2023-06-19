%%%% Continuous Gaussian Process-Based Search(GPS-C) Algorithm
clear,clc, close all;

% setting of the problem
fun_name = 'Fun_WeightedSphere';
 % get the information of the test problem
switch fun_name
        case 'Fun_WeightedSphere'
        num_vari=6; design_space=[-5.12,-5.12,-5.12,-5.12,-5.12,-5.12;5.12,5.12,5.12,5.12,5.12,5.12];   optimum=0;
        % num_vari=5; design_space=[-5.12,-5.12,-5.12,-5.12,-5.12;5.12,5.12,5.12,5.12,5.12];   optimum=0;
    otherwise
        error('objective function is not defined!')
end  
%--------------------------------------------------------------------------

d=4;       %dimension of the problem
MaxSteps=600;      %Maximum iteration of the algorithm
s=10;       %number of solutions sampled in each iteration

s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

% Step 0 - initializaiton
Opti_Record = zeros(MaxSteps,1);
Opti_Record_true = zeros(MaxSteps,1);
Opti_Record_Point = zeros(MaxSteps,d);

% initial design points using Latin hypercube sampling method
num_initial = 600;
Sample_ini = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
G_ini = feval(fun_name,Sample_ini);

% fit a new model
B = ones(size(Sample_ini,1)^1,1);       % basis function matrix at design points
skriging_model = SKfit(Sample_ini(:,1:d),G_ini,B,2);   %fit a kriging model

% randomly generate s points
Sample = 10.24*(rand(s,d)-0.5*ones(s,d));

% take one observation
SSample = [];
G = [];
G = feval(fun_name,Sample);
SSample = Sample;
SSample(:,d+1) = ones(size(SSample,1),1);
SSample(:,d+2) = G;
SSample(:,d+3) = skriging_model.sigma2 * ones(size(SSample,1),1);

% matrix inversion
sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d));
sigma_old = sigma;

% find current sample best
B = ones(size(SSample,1)^1,1);
SSample(:,d+4) = SKpredict(skriging_model,SSample(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
[g_star_sample_best,max_location] = max(SSample(:,d+4));
x_sample_best = SSample(max_location,1:d);

% find best solution
B = ones(size(1,1)^1,1);
[g_star_1,x_star_1] = find_current_best(skriging_model,SSample(:,1:d),B,SSample(:,d+2),x_sample_best,sigma);

% Step2 - iteration
step = 1;

Opti_Record(1) = g_star_1;
Opti_Record_Point(1,1:d) = x_star_1;

while step<MaxSteps
    % sample s points from fitted distribution
    tic;
    Sample = Sampling(skriging_model,SSample,g_star_1,x_star_1,s,sigma);
    toc;
    G = [];
    G = feval(fun_name,Sample);
    
    % Merge to SSample
    for i = 1:size(Sample,1)
            m = size(SSample,1);
            SSample(m+1,d+2) = G(i);
            SSample(m+1,1:d+1) = [Sample(i,:),1];
            SSample(m+1,d+3) = skriging_model.sigma2;
    end
    % matrix inversion
    tic;
    sigma = Matrix_inverse_1(skriging_model,SSample(:,d+3),SSample(:,1:d),sigma_old);
    sigma_old = sigma;
    toc;
    
    % find current sample best
    B = ones(size(SSample,1)^1,1);
    SSample(:,d+4) = SKpredict(skriging_model,SSample(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
    [g_star_sample_best,max_location] = max(SSample(:,d+4));
    x_sample_best = SSample(max_location,1:d);
    
    % find best solution
    B = ones(size(1,1)^1,1);
    [g_star_1,x_star_1] = find_current_best(skriging_model,SSample(:,1:d),B,SSample(:,d+2),x_sample_best,sigma);
    
    step = step+1;
    % fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f] %.2f \n',step,x_star_1(1),x_star_1(2),x_star_1(3),x_star_1(4),x_star_1(5),x_star_1(6),x_star_1(7),x_star_1(8),x_star_1(9),x_star_1(10),g_star_1);
    fprintf('%d  [%.2f, %.2f, %.2f, %.2f] %.2f \n',step,x_star_1(1),x_star_1(2),x_star_1(3),x_star_1(4),g_star_1);
    % fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f£¬ %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f] %.2f \n',step,x_star_1(1),x_star_1(2),x_star_1(3),x_star_1(4),x_star_1(5),x_star_1(6),x_star_1(7),x_star_1(8),x_star_1(9),x_star_1(10),x_star_1(11),x_star_1(12),x_star_1(13),x_star_1(14),x_star_1(15),x_star_1(16),x_star_1(17),x_star_1(18),x_star_1(19),x_star_1(20),g_star_1);
    % fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f] %.2f \n',step,x_star_1(1),x_star_1(2),x_star_1(3),x_star_1(4),x_star_1(5),g_star_1);
    Opti_Record(step) = g_star_1;
    Opti_Record_Point(step,1:d) = x_star_1;

end

for si=1:MaxSteps
    Opti_Record_true(si) = Fun_WeightedSphere_true(Opti_Record_Point(si,:));
end

