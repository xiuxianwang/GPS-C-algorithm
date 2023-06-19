% resample a point in the sampled point set
function p = Resampling(k_prime,X_prime,T)
    d=2;
    
    size_X = size(X_prime,1);
    T_prime = T/(log(k_prime+1));
    summation = sum(exp(X_prime(:,d+2)/T_prime));
    
    u = rand(1);
    
    cdf(1) = exp(X_prime(1,d+2)/T_prime)/summation;
    if u < cdf(1)
        p = 1;
    end
    for i=2:size_X
        cdf(i) = cdf(i-1) + exp(X_prime(i,d+2)/T_prime)/summation;
        if u < cdf(i)
            p = i;
            break
        end
    end

end