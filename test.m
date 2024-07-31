
clc

xdata = -4:0.1:5;
N = numel(xdata);

noise1 = 0.05*(2*rand(1,N)-1);

x=xdata;
yexp1 = 0.01*x.^2 + 0.1*x - 0.2 + noise1;



figure
plot(xdata, yexp1)


%%
clc

model_1 = @(a, b, c, x) a*x.^2 + b*x + c;


model = @(v) [model_1(v(1), v(2), v(3), xdata) - yexp1];


% out = model([0.2 1]);
% plot(xdata, out)



Lower = [ -inf     -inf   -inf];
Start = [  0.1        1    0];
Upper = [  inf      inf   inf];


options = optimoptions('lsqnonlin', ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations', 800000, ...
    'FunctionTolerance', 1E-9, ...
    'Algorithm','trust-region-reflective', ... %levenberg-marquardt trust-region-reflective
    'MaxIterations', 5000, ...
    'StepTolerance', 1e-9, ...
    'PlotFcn', '', ... %optimplotresnorm optimplotstepsize OR ''  (for none)
    'Display', 'iter', ... %final off iter
    'FiniteDifferenceStepSize', 1e-9, ...
    'CheckGradients', true, ...
    'DiffMaxChange', 0.01);

% [vestimated,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(ModelFunction, Start, Lower, Upper, options);
[vout,resnorm,residual,~,~,~,jacobian] = lsqnonlin(model, Start, Lower, Upper, options);



Hess = jacobian'*jacobian;

DataSize = size(residual,2);
ParSize = size(vout,2);

sigmSqr = (residual*residual')/(DataSize-ParSize);

VarCovar = inv(Hess)*sigmSqr;

for i=1:ParSize
    for j=1:ParSize
        Correl(i,j) = VarCovar(i,j)./(VarCovar(i,i)*VarCovar(j,j))^0.5;
    end
end
imagesc(Correl);


for k=1:ParSize
    error(k) = (full(VarCovar(k,k))).^0.5;
end

vout
error*3 %3 сигма














































