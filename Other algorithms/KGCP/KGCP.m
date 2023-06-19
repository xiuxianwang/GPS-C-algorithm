%--------------------------------------------------------------------------
clearvars; close all;
% setting of the problem
fun_name = 'Fun_Sixhump';
 % get the information of the test problem
switch fun_name
    case 'Fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum= 1.031628;
    case 'Fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= -0.397887;
    case 'Fun_Hills'
        num_vari=2; design_space=[0,0;100,100];                  optimum=20;
    case 'Fun_Rosenbrock'
        num_vari=10; design_space=[-10,-10,-10,-10,-10,-10,-10,-10,-10,-10;10,10,10,10,10,10,10,10,10,10];            optimum=0;
    otherwise
        error('objective function is not defined!')
end  
 %--------------------------------------------------------------------------

s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);
% Initialization
% the number of initial design points
num_initial =20;
% the number of total allowed design points
max_evaluation = 800;
%--------------------------------------------------------------------------
% generate a grid to evaluate KG
% % grid for sixhump prolem
grid = zeros(51*51,2);
val = zeros(51*51,2);
count = 1;
for i = 1:1:51
    for j = 1:1:51
        grid(count,:) = [-2+0.08*(i-1),-2+0.08*(j-1)];
        count = count +1;
    end
end
% % grid for branin problem
%     grid = zeros(51*51,2);
%     val = zeros(51*51,1);
%     count = 1;
%     for i = 1:1:51
%         for j = 1:1:51
%             grid(count,:) = [-5+0.3*(i-1),0+0.3*(j-1)];
%             count = count +1;
%         end
%     end
% % grid for GPS problem
% grid = zeros(51*51,2);
% val = zeros(51*51,2);
% count = 1;
% for i = 1:1:51
%     for j = 1:1:51
%         grid(count,:) = [0+2*(i-1),0+2*(j-1)];
%         count = count +1;
%     end
% end
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% record the f_min in each iteration
f_min = zeros(max_evaluation - num_initial + 1,1);
% the current best solution
f_min(1) = max(sample_y);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
count_iteration = 1;
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(1), optimum);
%--------------------------------------------------------------------------
% the iteration
while evaluation <  max_evaluation
    % build (or rebuild) the initial Kriging model
    tic;
    if evaluation <= 200
        % build (or rebuild) the initial Kriging model
        B =  ones(size(sample_x,1),1);
        kriging_model = SKfit(sample_x,sample_y,B,2);
        kriging_model_old = kriging_model;
    else
        if evaluation == 200 + 50 * count_iteration
            B =  ones(size(sample_x,1),1);
            kriging_model = SKfit(sample_x,sample_y,B,2);
            kriging_model_old = kriging_model;
            count_iteration = count_iteration + 1;
        else
            B =  ones(size(sample_x,1),1);
            kriging_model = SKfit_no_update(sample_x,sample_y,B,kriging_model_old,2);
        end
    end
    toc;
%--------------------------------------------------------------------------    
    % normalize sample_x first
    minX = kriging_model.minX;
    maxX = kriging_model.maxX;
    theta = kriging_model.theta;
    tau2 = kriging_model.tausquared;
    K1 = size(sample_x,1);     % number of prediction points
    sample_x_normalize = (sample_x - repmat(minX,K1,1)) ./ repmat(maxX-minX,K1,1);
    % calculate covariance matrix used in KGCP
    % covariance matrix
    K_sample_x = zeros(size(sample_x,1),size(sample_x,1));
    for i = 1:size(sample_x,1)
        x_sample_normalize = sample_x_normalize(i,:);
        K_sample_x(i,:) = cov_vector_speed(x_sample_normalize',sample_x_normalize',theta,tau2);
    end
    % covariance matrix plus noise variance matrix
    noise_variance = kriging_model.sigma2;
    K_sample_x_NoiseVar = K_sample_x + diag(noise_variance*ones(size(sample_x,1),1));
    % save the vector that depends only on V and y
    mu_0 = kriging_model.beta;
    K_Noise_Inv_y_mu = K_sample_x_NoiseVar \ (sample_y - mu_0*ones(size(sample_x,1),1));
    % KGCP value
    KGCP_val = @(x)KGCP_Compute_pro(x,sample_x_normalize',K_sample_x,K_sample_x_NoiseVar,K_Noise_Inv_y_mu,theta,tau2,mu_0,noise_variance);
    KGCP_val_2 = @(x)arrayfun(KGCP_val,x);
    tic;
    grid_normalize = (grid - repmat(minX,size(grid,1),1)) ./ repmat(maxX-minX,size(grid,1),1);
    ss = num2cell(grid_normalize,2); 
    val= feval(KGCP_val_2,ss);
    [value,index] = max(val);
    myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
    start = grid_normalize(index,:);
    lb = [max(start(1,1) - 0.02,grid_normalize(1,1)),max(start(1,2)-0.02,grid_normalize(1,2))];
    ub = [min(start(1,1) + 0.02,grid_normalize(2601,1)),min(start(1,2)+0.02,grid_normalize(2601,2))];
    best_x_normalize = fmincon(@(x)  KGCP_Compute_pro_fmincon(x,sample_x_normalize',K_sample_x,K_sample_x_NoiseVar,K_Noise_Inv_y_mu,theta,tau2,mu_0,noise_variance),...
                 start,[],[],[],[],lb,ub,[],myopt);
    best_x = best_x_normalize.* repmat(maxX-minX,1,1) + repmat(minX,1,1);
    toc;
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    evaluation = evaluation + 1;
    % add the new point to design set
    x_ind = find(sample_x(:,1) == best_x(1,1));
    y_ind = find(sample_x(:,2) == best_x(1,2));
    if (isempty(x_ind))
        sample_x = [sample_x;best_x];
        sample_y = [sample_y;best_y];
    elseif (isempty( find(sample_x(x_ind,2) == best_x(1,2) )))
        sample_x = [sample_x;best_x];
        sample_y = [sample_y;best_y];
    else
        ind_ind = find(sample_x(x_ind,2) == best_x(1,2));
        sample_y(x_ind(ind_ind,1),1) = 1/2*(sample_y(x_ind(ind_ind,1),1) + best_y);
    end
    B =  ones(size(sample_x,1),1);
    Yhat = SKpredict(kriging_model,sample_x,B);
    % update some parameters
    iteration = iteration + 1;
    [max_num,max_index] = max(Yhat);
    f_min(iteration+1) = Fun_Sixhump_free(sample_x(max_index,:));
     % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(iteration+1), optimum);
end





