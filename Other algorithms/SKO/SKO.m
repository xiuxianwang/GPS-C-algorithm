%--------------------------------------------------------------------------
clearvars; close all;
% setting of the problem
fun_name = 'Fun_Hills';
 % get the information of the test problem
switch fun_name
    case 'Fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=-1.031628;
    case 'Fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'Fun_Hills'
        num_vari=2; design_space=[0,0;100,100];                  optimum=-20;
    otherwise
        error('objective function is not defined!')
end  
 %--------------------------------------------------------------------------
s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);
% Initialization
% the number of initial design points
num_initial = 20;
% the number of total allowed design points
max_evaluation = 800;
% % grid for branin problem
% grid = zeros(51*51,2);
% val = zeros(51*51,1);
% count = 1;
% for i = 1:1:51
%     for j = 1:1:51
%         grid(count,:) = [-5+0.3*(i-1),0+0.3*(j-1)];
%         count = count +1;
%     end
% end
 % grid for sixhump prolem
%     grid = zeros(51*51,2);
%     val = zeros(51*51,1);
%     count = 1;
%     for i = 1:1:51
%         for j = 1:1:51
%             grid(count,:) = [-2+0.08*(i-1),-2+0.08*(j-1)];
%             count = count +1;
%         end
%     end
% % grid for Hills problem
grid = zeros(51*51,2);
val = zeros(51*51,1);
count = 1;
for i = 1:1:51
    for j = 1:1:51
        grid(count,:) = [0+2*(i-1),0+2*(j-1)];
        count = count +1;
    end
end
%--------------------------------------------------------------------------
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% record the f_min in each iteration
f_min = zeros(max_evaluation - num_initial + 1,1);
% the current best solution
f_min(1) = min(sample_y);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
count_iteration = 1;
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(1), optimum);
%--------------------------------------------------------------------------
% the iteration
while evaluation <  max_evaluation
    if evaluation <= 200
        % build (or rebuild) the initial Kriging model
        B =  ones(size(sample_x,1),1);
        kriging_model = SKfit(sample_x,sample_y,B,2);
        kriging_model_old = kriging_model;
    else
        B =  ones(size(sample_x,1),1);
        kriging_model = SKfit_no_update(sample_x,sample_y,B,kriging_model_old,2);
    end
    toc;
    % find the x^** using u(x)
    [Yhat,MSE] = SKpredict(kriging_model,sample_x,B);
    [max_num,max_index] = max(-Yhat - MSE.^(1/2));
    Yhat_star= Yhat(max_index);
    % the Expected Improvement criterion
    infill_criterion = @(x)Infill_Standard_EI(x,kriging_model,Yhat_star);
    % find the highest EI value using a grid evaluation + fmincon
    tic;
    val = feval(infill_criterion,grid);
    [value,index] = min(val);
    myopt = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',50);
    start = grid(index,:);
    lb = [max(start(1,1) - 2,0),max(start(1,2)-2,0)];
    ub = [min(start(1,1) + 2,100),min(start(1,2)+2,100)];
    best_x = fmincon(@(x)  Infill_Standard_EI_fmincon(x,kriging_model,Yhat_star),...
             start,[],[],[],[],lb,ub,[],myopt);
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
        count_int = 0;
    elseif (isempty( find(sample_x(x_ind,2) == best_x(1,2) )))
        sample_x = [sample_x;best_x];
        sample_y = [sample_y;best_y];
        count_int = 0;
    else
        ind_ind = find(sample_x(x_ind,2) == best_x(1,2));
        sample_y(x_ind(ind_ind,1),1) = 1/2*(sample_y(x_ind(ind_ind,1),1) + best_y);
        count_int = 1;
    end
    if (count_int ==0)
        B =  ones(size(best_x,1),1);
        [sample_y_est,MSE_est] = SKpredict(kriging_model,best_x,B);
        Yhat = [Yhat;sample_y_est];
    end
    % update some parameters
    iteration = iteration + 1;
    [min_num,min_index] = min(Yhat);
    f_min(iteration+1) = Fun_Hills_free(sample_x(min_index,:));
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(iteration+1), optimum);
end

    



