function obj = Fun_Hills_free(x)
%--------------------------------------------------------
% The Hills problem without noise
% range[0,100]x[0,100], global = [90,90]
% fmax = 20
%--------------------------------------------------------
obj = 10*(sin(0.05*pi.*x(:,1))).^6./(2.^(2*((x(:,1)-90)/80).^2))+...
   10*(sin(0.05*pi.*x(:,2))).^6./(2.^(2*((x(:,2)-90)/80).^2));

end













