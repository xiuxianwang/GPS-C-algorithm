% many local
function y = Fun_Rosenbrock(x)

y=zeros(size(x,1),1);
% The Rosenbrock example
% range[-10,10]^10, global = [1,1,...,1]
for i=1:9
    y = y+(1-x(:,i)).^2+100.*(x(:,i+1)-x(:,i).^2).^2;
end
y = -10^(-4)*y;
   


% add noise N(0,0.1)
y = y + normrnd(zeros(size(x,1),1),min(0.1*(ones(size(x,1),1)+abs(y)),5*ones(size(x,1),1)),[size(x,1),1]);

end











