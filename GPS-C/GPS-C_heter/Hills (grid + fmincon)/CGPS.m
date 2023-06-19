%%%% Continuous Gaussian Process-Based Search(GPS-C) Algorithm
clear,clc, close all;

% setting of the problem
fun_name = 'Fun_Hills';
 % get the information of the test problem
switch fun_name
    case 'Fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=-1.031628;
    case 'Fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'Fun_Hills'
        num_vari=2; design_space=[0,0;100,100];                optimum=20;
    case 'Fun_Rosenbrock'
        num_vari=10; design_space=[-10,-10,-10,-10,-10,-10,-10,-10,-10,-10;10,10,10,10,10,10,10,10,10,10];            optimum=0;
    otherwise
        error('objective function is not defined!')
end  
%--------------------------------------------------------------------------


d=2;       %dimension of the problem
MaxSteps=190;      %Maximum iteration of the algorithm
r0 = 2;     % initial bandwidth of the kernel
beta = 0.009;    % contracting rate of the kernel bandwidth
s=10;       %number of solutions sampled in each iteration

% generate uniform grids in the feasible region for evaluation
count = 0;
A = zeros(51*51,2);
for i=1:51
    for j = 1:51
        count = count+1;
        A(count,1) = (i-1)*2;
        A(count,2) = (j-1)*2;
    end
end

%% The main steps of the GPS-C algorithm

s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);
% Step 0 - initializaiton

Opti_Record = zeros(MaxSteps,1);

% initial design points using Latin hypercube sampling method
num_initial = 100;
Sample_ini = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
G_ini = feval(fun_name,Sample_ini);

% randomly sample s points in the feasible region (optional)
Sample = 100*rand(s,2);

% take and record the observations
SSample = [];
dist_matrix = [];
G = [];
G = feval(fun_name,Sample);     % take observations of the randomly sampled points
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
    if (lambda_hat<0.2)    % set a lower bound for the variance of simulation noise
        SSample(i,d+3) = 0.2;
    else
         SSample(i,d+3) = lambda_hat;
    end
end

% fit a krging model
B = ones(size(SSample,1)^1,1);       % basis function matrix at design points
skriging_model = SKfit2(SSample(:,1:d),SSample(:,d+2),B, SSample(:,d+3), 2);   %fit a kriging model (this kriging method differs from the Function SKfit, since it intensionally returns larger tau2hat for better exploration in early stage of the algorithm)

% matrix inverse for further calculating the conditional mean and variance
sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d));

% find current best solution using grid search + fmincon
B = ones(size(A,1)^1,1);
[Y_hat,sigma]  = SKpredict(skriging_model,A(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
[start_value,max_location] = max(Y_hat);    % find the optimal solution of the grid
start = A(max_location,1:d);    % use the grid-best solution as the start point
myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
lb = [max(start(1,1) - 2,0),max(start(1,2)-2,0)];
ub = [min(start(1,1) + 2,100),min(start(1,2)+2,100)];
B = ones(1^1,1);
x_star = fmincon(@(x) SKpredict_fmincon(skriging_model,x,SSample(:,1:d),B,SSample(:,d+2),sigma),...
                     start,[],[],[],[],lb,ub,[],myopt);
g_star_estimate = SKpredict(skriging_model,x_star,SSample(:,1:d),1,SSample(:,d+2),sigma);
g_star = Fun_Hills_free(x_star);    % calulate the true objective function value

% Step2 - iteration
step = 1;
Opti_Record(1) = g_star;

while step<MaxSteps
    % sample s points from fitted distribution
    tic;
    Sample = Sampling(skriging_model,SSample,g_star_estimate,x_star,s,sigma);
    toc;
    % take one observation for each newly sampled point
    G = [];
    G = feval(fun_name,Sample);
    rk =r0 * (k + 1)^(-beta);  % update the bandwidth of the kernel
    
    % Merge to SSample
    for i = 1:size(Sample,1)
            m = size(SSample,1);
            SSample(m+1,d+2) = G(i);
            SSample(m+1,1:d+1) = [Sample(i,:),1]; 
    end
    % update the distance matrix
    k = size(SSample,1);
    tic;
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
        if (lambda_hat<0.2)
            SSample(i,d+3) = 0.2;
        else
            SSample(i,d+3) = lambda_hat;
        end
    end
    toc;
    % continue updating the Gaussian process parameters in early stage
    if (step==20)
        B = ones(size(SSample,1)^1,1);       % basis function matrix at design points
        skriging_model = SKfit2(SSample(:,1:d),SSample(:,d+2),B, SSample(:,d+3),2);   %fit a kriging model
    end
    if (step==40)
        B = ones(size(SSample,1)^1,1);       % basis function matrix at design points
        skriging_model = SKfit(SSample(:,1:d),SSample(:,d+2),B, SSample(:,d+3),2);   %fit a kriging model
    end
    % matrix inversion
    tic;
    sigma = Matrix_inverse(skriging_model,SSample(:,d+3),SSample(:,1:d));
    toc;
    
    % find current best solution using grid search + fmincon
    B = ones(size(A,1)^1,1);
    [Y_hat,sigma]  = SKpredict(skriging_model,A(:,1:d),SSample(:,1:d),B,SSample(:,d+2),sigma);
    [start_value,max_location] = max(Y_hat);    % find the optimal solution of the grid
    start = A(max_location,1:d);    % use the grid-best solution as the start point
    myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
    lb = [max(start(1,1) - 2,0),max(start(1,2)-2,0)];
    ub = [min(start(1,1) + 2,100),min(start(1,2)+2,100)];
    B = ones(1^1,1);
    x_star = fmincon(@(x) SKpredict_fmincon(skriging_model,x,SSample(:,1:d),B,SSample(:,d+2),sigma),...
                         start,[],[],[],[],lb,ub,[],myopt);
    g_star_estimate = SKpredict(skriging_model,x_star,SSample(:,1:d),1,SSample(:,d+2),sigma);
    g_star = Fun_Hills_free(x_star);    % calulate the true objective function value
    
    step = step+1;    % update the iteration count
    fprintf('%d  [%.2f, %.2f] %.2f \n',step,x_star(1),x_star(2),g_star);    % output the current best solution
    
    Opti_Record(step) = g_star;     % record the objective value of current best solution

end



% %%%%% plot the evaluated points
% figure, hold on;
% axis equal;
% for i=1:100*s
%     scatter(SSample(i,1),SSample(i,2),70,'r.');
% end
% xlim([0,100]);
% ylim([0,100]);
% ylabel('GPS-C','Fontsize',18);
% title('1000','Fontsize',18);
% set(gca,'Ytick',0:20:100);
% box on;



