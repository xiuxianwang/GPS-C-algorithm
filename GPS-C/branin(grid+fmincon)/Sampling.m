%%%%% construct sampling distribution and sample solution
function Sample = Sampling(skriging_model_2,SSample,g_star,x_star,s,sigma)

d=2;       %dimension
mean_LB = -10;       %the lower bound of the conditional mean
var_LB = 0.01;      %the lower bound of the conditional variance

%%% Markov Chain Coordinate Sampling (MCCS)
T = 100;
Sample = zeros(s,d);
for i=1:s   %number of solutions sampled
    y = x_star;
    for t = 1:T
        % Sample Uniformly an interger I from 1 to d
        I = unidrnd(d);
        % Sample discrete j uniformly from [0,y(I)-0.1]u[y(I)+0.1,10]
        j = y(I);
        while j == y(I)
            new_cand = 15*(rand(1,2)-repmat([1/3,0],1,1));  
            j = new_cand(I);
        end
        z = y;
        z(I) = j;
        % use z or not?
               
        x = y;
        B = ones(size(x,1)^1,1);       % basis function matrix at design points
        
        E_Y = SKpredict(skriging_model_2,x,SSample(:,1:d),B,SSample(:,d+2),sigma);
        Var_Y = CalculateMSE(skriging_model_2,x,sigma,SSample(:,1:d),B);
        
        if(E_Y<mean_LB)
            E_Y=mean_LB;
        end
        if(Var_Y<var_LB)
            Var_Y=var_LB;
        end
        
        P_Y_y = (1 - normcdf(g_star,E_Y,Var_Y^0.5)); % returns the normal cdf at x for the mean and standard deviation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = z;       
        B = ones(size(x,1)^1,1);       % basis function matrix at design points
        E_Y = SKpredict(skriging_model_2,x,SSample(:,1:d),B,SSample(:,d+2),sigma);
        Var_Y = CalculateMSE(skriging_model_2,x,sigma,SSample(:,1:d),B);  
        
        if(E_Y<mean_LB)
            E_Y=mean_LB;
        end
        if(Var_Y<var_LB)
            Var_Y=var_LB;
        end
        
        P_Y_z = (1 - normcdf(g_star,E_Y,Var_Y^0.5)); % returns the normal cdf at x for the mean and standard deviation
        
        if rand <= P_Y_z/P_Y_y
            y = z;
        end     
    end
    Sample(i,:) = y;
end




