function y = Fun_Rosenbrock(x,n)


y_temp = 0;
y_sum = 0;
for i=1:n
% The Rosenbrock example
% range[-10,10]^10, global = [1,1,...,1]
    y=0;    
    for i=1:9
        y = y+(1-x(i))^2+100*(x(i+1)-x(i)^2)^2;
    end
    y_temp = -10^(-4)*y;
   
% add noise N(0,0.1)
    y_temp = y_temp+normrnd(0,min(0.1*(1+abs(y_temp)),5));
    y_sum = y_sum +y_temp;
end
y = y_sum/n;
end














