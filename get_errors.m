function errors = get_errors(vout, residual, jacobian)

ParSize = size(vout, 2);

VarCovar = full(inv(jacobian.'*jacobian)*var(residual));

for i=1:ParSize
    for j=1:ParSize
        Correl(i,j) = VarCovar(i, j)./(VarCovar(i, i)*VarCovar(j, j))^0.5;
    end
end
for i = 1:ParSize
    Correl(i,i) = 0;
end
% figure
% imagesc(Correl);
% load('correlation_colormap.mat');
% colormap(Correl_colormap);


for k=1:ParSize
    error(k) = (full(VarCovar(k,k))).^0.5;
end

% vout
errors = error*1; %1 сигма

end