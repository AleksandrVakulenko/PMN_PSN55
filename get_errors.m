function errors = get_errors(vout, residual, jacobian)

ParSize = size(vout, 2);

VarCovar = full(inv(jacobian.'*jacobian)*var(residual));

for i=1:ParSize
    for j=1:ParSize
        Correl(i,j) = VarCovar(i, j)./(VarCovar(i, i)*VarCovar(j, j))^0.5;
    end
end
% imagesc(Correl);


for k=1:ParSize
    error(k) = (full(VarCovar(k,k))).^0.5;
end

% vout
errors = error*3; %3 сигма

end