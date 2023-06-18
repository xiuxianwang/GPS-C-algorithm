% many local
function y = Fun_Rosenbrock_free(x)


y=0;
% The Rosenbrock example
% range[-10,10]^10, global = [1,1,...,1]
for i=1:9
    y = y+(1-x(:,i)).^2+100.*(x(:,i+1)-x(:,i).^2).^2;
end
y = -10^(-4)*y;
    
end

