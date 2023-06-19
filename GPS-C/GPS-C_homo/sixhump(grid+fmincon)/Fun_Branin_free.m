function obj = Fun_Branin_free(x)
%------------------------------------------------
% The Branin problem without noise
% range[-5,10]x[0,15], global = [9.42477796,2.47499998], [-3.14159265,12.27500000], [3.14159265,2.27500000]
% fmax = -0.397887357729739
%--------------------------------------------------------
a=1;
b=5.1/(4*pi*pi);
c=5/pi;
d=6;
h=10;
ff=1/(8*pi);
  x1 = x(:,1);
  x2 = x(:,2);
obj=a.*(x2-b.*x1.^2+c.*x1-d).^2+h.*(1-ff).*cos(x1)+h;
obj = -obj; % The Branin problem is a minimization problem
end
