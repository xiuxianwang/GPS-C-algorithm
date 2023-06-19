%%%%% construct sampling distribution and sample solution

function Sample = Sampling(x_star)

d = 2;

% sampling scheme 
v  = 2*(rand(1,d)-0.5*ones(1,d)); % arbitrarily choose a direction

% randomly choose one point along direction v
ind = find(v>0);
v_1 = v(ind);
x_star_1 = x_star(ind);
k = size(v_1,2);
upper = min((100*ones(1,k)-x_star_1)./v_1);

ind2 = find(v<0);
v_2 = v(ind2);
x_star_2 = x_star(ind2);
k = size(v_2,2);
lower = min((0*ones(1,k)-x_star_2)./v_2);

if k==0
    bound = upper;
elseif k==2
    bound = lower;
else
    bound = min(upper,lower);
end

s = bound*rand;

Sample = x_star + s*v;

end




