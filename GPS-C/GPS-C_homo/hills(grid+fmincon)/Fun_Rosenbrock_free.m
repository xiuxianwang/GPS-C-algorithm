function obj = Fun_Rosenbrock_free(x)
%--------------------------------------------------------
% The Rosenbrock problem without simulation noise
% range[-10,10]^10, global = [1,1,...,1]
% fmax = 0
%--------------------------------------------------------

obj=zeros(size(x,1),1);
for i=1:9
    obj = obj+(1-x(:,i)).^2+100.*(x(:,i+1)-x(:,i).^2).^2;
end

obj = -10^(-6)*obj;
obj = - obj;

end











