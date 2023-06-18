function y = fun_many_loc(x,skriging_model)

B = ones(size(x,1)^1,1);
y = SKpredict_sample_path(skriging_model,x,B);

% add noise N(0,0.5)
for i = 1:size(y,1)
    y(i) = y(i) + normrnd(0,0.5);
end

end













