%% intro


ColeCole=@(DEps, f, a, X) DEps./(1+(1i*2*pi.*X./f).^a);

ColeCole_4=@(DEps1, f1, a1, ...
             DEps2, f2, a2, ...
             DEps3, f3, a3, ...
             DEps4, f4, a4, X) ...
    ColeCole(DEps1, f1, a1, X) + ...
    ColeCole(DEps2, f2, a2, X) + ...
    ColeCole(DEps3, f3, a3, X) + ...
    ColeCole(DEps4, f4, a4, X);

Freal=@(v, X) real(ColeCole_4(v(1), v(2), v(3), ...
                                   v(4), v(5), v(6), ...
                                   v(7), v(8), v(9), ...
                                   v(10), v(11), v(12), ...
                                   X));


Fimag = @(v, X) -imag(ColeCole_4(v(1), v(2), v(3), ...
                                    v(4), v(5), v(6), ...
                                    v(7), v(8), v(9), ...
                                    v(10), v(11), v(12), ...
                                    X));


load('Data/coeff.mat')
load('Data/Temp.mat')



%% main fit loop

clc

Coeff_array = [];
Errors_array = [];

N = 20


filename = ['./Data/freq/T_' num2str(N) '.TXT'];
[freq, eps1, eps2] = importfile(filename);



relative_error = 0.6/100;
eps1_abs_error = eps1*relative_error;
eps2_abs_error = eps2*relative_error;

ModelFunction = @(v) [(Freal(v, freq) - eps1)./eps1_abs_error; ...
                      (Fimag(v, freq) - eps2)./eps2_abs_error]';

% ModelFunction = @(v) [Fimag(v, xdata) - ydata(:,2)]';



Start = coeff(N, :);

% dEps
Lower([1,4,7,10]) = 0;
Upper([1,4,7,10]) = 50000;

%freq
Lower([2,5,8,11]) = coeff(N, [2,5,8,11])*0.8;
Upper([2,5,8,11]) = coeff(N, [2,5,8,11])*1.2;

% alpha
Lower([3,6,9,12]) = [0, 1, 0, 0];
Upper([3,6,9,12]) = [1, 1, 1, 1];


options = optimoptions('lsqnonlin', ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations', 800000, ...
    'FunctionTolerance', 1E-12, ...
    'Algorithm','trust-region-reflective', ... %levenberg-marquardt trust-region-reflective
    'MaxIterations', 5000, ...
    'StepTolerance', 1e-12, ...
    'PlotFcn', '', ... %optimplotresnorm optimplotstepsize OR ''  (for none)
    'Display', 'final', ... %final off iter
    'FiniteDifferenceStepSize', 1e-12, ...
    'CheckGradients', true, ...
    'DiffMaxChange', 0.001, ...
    'OptimalityTolerance', 1e-12);

[vout, resnorm, residual, ~, ~, ~, jacobian] = lsqnonlin(ModelFunction, Start, Lower, Upper, options);



errors = get_errors(vout, residual, jacobian);

Coeff_array(N, :) = vout;
Errors_array(N, :) = errors;

PRC_1 = errors./vout*100;
%% look at single fit (last used)


eps1_model = Freal(vout, freq);
eps2_model = Fimag(vout, freq);

figure
hold on
% plot(freq, eps2, '-xb')
errorbar(freq, eps2, eps2*relative_error, '-xb')
plot(freq, eps2_model, '-r')
set(gca, 'xscale', 'log')

figure
hold on
% plot(freq, eps1, '-xb')
errorbar(freq, eps1, eps1*relative_error, '-xb')
plot(freq, eps1_model, '-r')
set(gca, 'xscale', 'log')

%%

resnum = numel(residual);
varnum = 12;
try_errors = 10.^linspace(-3, 4, 100);
chi2 = sum(residual'.^2./(try_errors).^2)/(resnum-varnum);

figure
plot(try_errors, chi2)
yline(1)
set(gca,'yscale','log')
set(gca,'xscale','log')





chi2 = sum(residual'.^2./(1).^2)/(resnum-varnum)


%%



plot(residual, '-x')
yline(0)


%%


h = histogram(residual, 20, 'Normalization', 'pdf');
values = h.Values;
bin_edges = h.BinEdges;
bins = (bin_edges(1:end-1)+bin_edges(2:end))/2;

hold on
plot(bins, values, 'x');


























