%%%% Single-replication based random search algorithm
clear,clc, close all;

%% parameter setting
d=10;
MaxSteps = 4000;
b = 1.1;
c = 0.5;
delta = 0.01;
L = 10;
T = 1;

%% Step 0-3
% Step 0 - initializaiton
s_rand = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s_rand);

i = 1;
k = 0;
X = [];

while k<MaxSteps
    
    k = k + 1;
    if  k == floor(i^b)
        if i ==1
            x = 20*(rand(1,10)-0.5*ones(1,10));
            X(1,1:d) = x;
            X(1,d+1) = L;
            X(1,d+2) = Fun_Rosenbrock(x,L);
            x_star = x;
            g_star = X(1,d+2);
            Opti_Record(1) = g_star;
            Opti_Record_Point(1,1:d) = x_star;
        else
            x = Sampling(x_star);
            f_L = Fun_Rosenbrock(x,L);
            if f_L >= Opti_Record(k-1) - delta
                count = size(X,1);
                X(count+1,1:d) = x;
                X(count+1,d+1) = L;
                X(count+1,d+2) = Fun_Rosenbrock(x,L);
            end
            for j = 1:size(X,1)
                if X(j,d+1) < ceil(i^c)
                    N_E = ceil(i^c) - X(j,d+1);
                    y_E = Fun_Rosenbrock(X(j,1:d),N_E);
                    X(j,d+2) = (X(j,d+1)*X(j,d+2)+N_E*y_E)/ceil(i^c);
                    X(j,d+1) = ceil(i^c);
                end
            end
        end
        i = i+1;
        k_prime = k;
        X_prime = X;
    else
        count = size(X,1);
        p = Resampling(k_prime,X_prime,T);
        y_E = Fun_Rosenbrock(X(p,1:d),1);
        X(p,d+2) = (X(p,d+1)*X(p,d+2)+y_E)/(X(p,d+1)+1);
        X(p,d+1) = X(p,d+1) + 1;       
    end
    
    [g_star,max_location] = max(X(:,d+2));
    x_star = X(max_location,1:d);           
    
    Opti_Record(k) = g_star;
    Opti_Record_Point(k,1:d) = x_star;
    Opti_Record_true(k) = Fun_Rosenbrock_free(x_star);
    
    fprintf('%d  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f] %.6f \n',k,x_star(1),x_star(2),x_star(3),x_star(4),x_star(5),x_star(6),x_star(7),x_star(8),x_star(9),x_star(10),g_star);
    
end


