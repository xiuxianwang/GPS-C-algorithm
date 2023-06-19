%%%% ASR algorithm
clear,clc, close all;

%% parameter setting
d=2;
MaxSteps = 2000;
b = 1.1;
c = 0.5;
delta = 1;
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
            x = 100*rand(1,2);
            X(1,1:d) = x;
            X(1,d+1) = L;
            X(1,d+2) = Fun_Hills(x,L);
            x_star = x;
            g_star = X(1,d+2);
            Opti_Record(1) = g_star;
            Opti_Record_Point(1,1:d) = x_star;
        else
            x = Sampling(x_star);
            f_L = Fun_Hills(x,L);
            if f_L >= Opti_Record(k-1) - delta
                count = size(X,1);
                X(count+1,1:d) = x;
                X(count+1,d+1) = L;
                X(count+1,d+2) = Fun_Hills(x,L);
            end
            for j = 1:size(X,1)
                if X(j,d+1) < ceil(i^c)
                    N_E = ceil(i^c) - X(j,d+1);
                    y_E = Fun_Hills(X(j,1:d),N_E);
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
        y_E = Fun_Hills(X(p,1:d),1);
        X(p,d+2) = (X(p,d+1)*X(p,d+2)+y_E)/(X(p,d+1)+1);
        X(p,d+1) = X(p,d+1) + 1;       
    end
    
    [g_star,max_location] = max(X(:,d+2));
    x_star = X(max_location,1:d);           
    
    Opti_Record(k) = g_star;
    Opti_Record_Point(k,1:d) = x_star;
    Opti_Record_true(k) = Fun_Hills_free(x_star);
    
    fprintf('%d  [%.2f, %.2f] %.2f \n',k,x_star(1),x_star(2),g_star);
    
end


% %%%%% plot the evaluated points
% figure, hold on;
% axis equal;
% box on;
% num = size(X,1);
% for i=1:num
%     plot(X(i,1),X(i,2),'r.');
% end
% axis([0 100 0 100])
% set(gca,'xtick',[0:20:100]);
% set(gca,'ytick',[0:20:100]);
% ylabel('ASR','Fontsize',18);
