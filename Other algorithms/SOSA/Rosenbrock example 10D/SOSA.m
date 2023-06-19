%%%% Single-replication based random search algorithm
clear,clc, close all;

%% parameter setting
d=10;
MaxSteps = 4000;
gamma = 0.91;
beta = 0.009;
s = 0.9;
kappa = 1;

%% Step 0 - initializaiton
s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);
Sk = [];
Gk = [];
dist_matrix = [];

% randomly sample the first point
Sk(1,1:d) = 20*(rand(1,10)-0.5*ones(1,10));
Gk(1) = Fun_Rosenbrock(Sk(1,1:d));
l(1) = 1;
dist_matrix(1,1) = 0;

x_star = Sk(1,1:d);
g_star = Gk(1);
f_hat(1) = Gk(1);
x_star_1 = Sk(1,1:d);

Opti_Record(1) = g_star;
Opti_Record_Point(1,1:d) = x_star;
Opti_Record_true(1) = Fun_Rosenbrock_free(x_star);

step = 1;

%% Step 0 - generate new design point

while step<MaxSteps
    
% sample a new points 
Sample = Sampling(x_star_1);
Sk(step+1,1:d) = Sample;
Gk(step+1) = Fun_Rosenbrock(Sample);
rk = (step + 1)^(-beta);

k = size(Sk,1);

for i = 1:k
    dist_matrix(i,k) = norm(Sk(i,1:d)-Sk(k,1:d));
    dist_matrix(k,i) = norm(Sk(i,1:d)-Sk(k,1:d));
end

for i = 1:k-1
    if dist_matrix(i,k)<rk
        f_hat(i) = (f_hat(i)*l(i) + Gk(k))/(l(i)+1);
        l(i) = l(i)+1;
    end
end

ind = [];
ind = find(dist_matrix(k,:)<rk);
l(k) = size(ind,2);
f_hat(k) = sum(Gk(ind))/l(k);

i_n = floor((step+1)^(0.9));

[g_star,max_location] = max(f_hat(1:i_n));
x_star = Sk(max_location,1:d);

[g_star_1,max_location] = max(f_hat(1:k));
x_star_1 = Sk(max_location,1:d);

step = step+1;

Opti_Record(step) = g_star;
Opti_Record_Point(step,1:d) = x_star;
Opti_Record_true(step) = Fun_Rosenbrock_free(x_star);

fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f] %.2f \n',step,x_star(1),x_star(2),x_star(3),x_star(4),x_star(5),x_star(6),x_star(7),x_star(8),x_star(9),x_star(10),g_star);

end
