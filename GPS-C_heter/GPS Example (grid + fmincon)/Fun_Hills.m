% many local
function obj = Fun_Hills(x)

% The Hills example
% range[0,100]x[0,100], global = [90,90]
obj = 10*(sin(0.05*pi.*x(:,1))).^6./(2.^(2*((x(:,1)-90)/80).^2))+...
   10*(sin(0.05*pi.*x(:,2))).^6./(2.^(2*((x(:,2)-90)/80).^2));

obj = obj + normrnd(zeros(size(x,1),1),0.5*sqrt(abs(obj)),[size(x,1),1]);

end










