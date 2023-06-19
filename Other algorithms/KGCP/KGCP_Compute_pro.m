function KGCP_val = KGCP_Compute_pro(xss,V,K_sample_x,K_sample_x_NoiseVar,K_Noise_Inv_y_mu,theta,tau2,mu0,noise_variance)


x = xss{1}';
k_V_V = K_sample_x;
k_V_V_w_NoiseVar = K_sample_x_NoiseVar;
V_Noise_Inv_y_mu = K_Noise_Inv_y_mu;

n = size(V,2);
ab = zeros(n+1,2);

k_v_V = cov_vector_speed(x,V,theta,tau2);
kn_v_v = tau2 - k_v_V * (k_V_V_w_NoiseVar \ k_v_V');

ab(1:n,1) = repmat(mu0,n,1) + k_V_V * V_Noise_Inv_y_mu;
kn_v_V = k_v_V' - k_V_V * (k_V_V_w_NoiseVar \ k_v_V');
ab(1:n,2) = kn_v_V / sqrt(kn_v_v + noise_variance);

ab(n+1,1) = mu0 + k_v_V * V_Noise_Inv_y_mu;
ab(n+1,2) = kn_v_v / sqrt(kn_v_v + noise_variance);

ab = sortrows(ab,2);

c = inf(1,n+2);
c(1) = -inf;
A = 0;
for i = 1:n
    if ab(i+1,2) == ab(i,2) && ab(i+1,1) <= ab(i,1)
        continue;
    else
        loopdone = 0;
        while loopdone == 0
            j = A(end);
            c(j+1+1) = (ab(j+1,1) - ab(i+1,1)) / (ab(i+1,2) - ab(j+1,2));
            if length(A) ~= 1 && c(j+1+1) <= c(A(end-1)+1+1)
                A(end) = [];
            else
                A = [A,i];
                loopdone = 1;
            end
        end
    end
end

a_new = zeros(1,length(A));
b_new = a_new;
c_new = inf(1,length(A)+1);
c_new(1) = -inf;
for i = 1:length(A)
    a_new(i) = ab(A(i)+1,1);
    b_new(i) = ab(A(i)+1,2);
    c_new(i+1) = c(A(i)+1+1);
end

KGCP_val = a_new .* (normcdf(c_new(2:end)) - normcdf(c_new(1:end-1))) + ...
            b_new .* (normpdf(c_new(1:end-1)) - normpdf(c_new(2:end)));
KGCP_val = sum(KGCP_val);

KGCP_val =  KGCP_val - max(ab(:,1));