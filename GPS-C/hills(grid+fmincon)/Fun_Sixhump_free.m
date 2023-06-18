function obj = Fun_Sixhump_free(x)
%--------------------------------------------------------
% 6-Hump Camel problem without noise
% range[-2,2]x[-2,1], global = [0.08984201,-0.71265640], [-0.08984201,0.71265640]           
% fmax = 1.0316284535
%---------------------------------------------------------
x1 = x(:,1);
x2 = x(:,2);
obj =(4-2.1.*x1.^2+x1.^4./3).*x1.^2+x1.*x2+(-4+4.*x2.^2).*x2.^2; 
obj = -obj;
end
