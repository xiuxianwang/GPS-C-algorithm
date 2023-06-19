function obj = Fun_Branin_free(x)
%------------------------------------------------
% Branin Test Function for Nonlinear Optimization
% Taken from "Towards Global Optimisation 2",edited by L.C.W. Dixon and G.P.
% Szego, North-Holland Publishing Company, 1978. ISBN 0 444 85171 2
%
% -5 <= x1 <= 10                  
%  0 <= x2 <= 15                  
% fmin = 0.397887357729739
% xmin =   9.42477796   -3.14159265  3.14159265
%          2.47499998   12.27500000  2.27500000
%---------------------------------------------------------
a=1;
b=5.1/(4*pi*pi);
c=5/pi;
d=6;
h=10;
ff=1/(8*pi);
  x1 = x(:,1);
  x2 = x(:,2);
obj=a.*(x2-b.*x1.^2+c.*x1-d).^2+h.*(1-ff).*cos(x1)+h;

obj = -obj;
