% notice that ASR can take multiple observations at each sampled point
function y = Fun_Hills(x,n)

y_temp = 0;
y_sum = 0;
for i=1:n

% The GPS algorithm example
% range[0,100]x[0,100], global = [90,90]
y_temp = 10*(sin(0.05*pi*x(1)))^6/(2^(2*((x(1)-90)/80)^2))+...
   10*(sin(0.05*pi*x(2)))^6/(2^(2*((x(2)-90)/80)^2));


% add noise
 y_temp = y_temp+normrnd(zeros(size(x,1),1),0.5*sqrt(abs(y_temp)),[size(x,1),1]);
 y_sum = y_sum +y_temp;
end
y = y_sum/n;  %return the mean value
 
end







