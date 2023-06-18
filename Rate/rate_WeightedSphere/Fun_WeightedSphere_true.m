% true objective function value
function y = Fun_WeightedSphere_true(x)


y=0;
% The Weighted sphere example
for i=1:4  % the dimension of this problem can be changed
    y = y+i*x(:,i).^2;
end   
y = -y;

end










