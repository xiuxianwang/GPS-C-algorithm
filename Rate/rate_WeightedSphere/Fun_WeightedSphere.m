% many local
function y = Fun_WeightedSphere(x)

y=zeros(size(x,1),1);
% The Weighted Sphere example
for i=1:4  % the dimension of this problem can be changed
    y = y+i*x(:,i).^2;
end   
y = -y;

% add noise N(0,0.1)
 y = y + normrnd(zeros(size(x,1),1),0.1*(ones(size(x,1),1)),[size(x,1),1]);

end











