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

for N = 2:49
    disp(N)

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
Lower([2,5,8,11]) = coeff(N, [2,5,8,11])*0.5;
Upper([2,5,8,11]) = coeff(N, [2,5,8,11])*2.0;

% alpha
Lower([3,6,9,12]) = [0, 1, 0, 0];
Upper([3,6,9,12]) = [1, 1, 1, 1];


options = optimoptions('lsqnonlin', ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations', 80000, ...
    'FunctionTolerance', 1E-9, ...
    'Algorithm','trust-region-reflective', ... %levenberg-marquardt trust-region-reflective
    'MaxIterations', 5000, ...
    'StepTolerance', 1e-10, ...
    'PlotFcn', '', ... %optimplotresnorm optimplotstepsize OR ''  (for none)
    'Display', 'final', ... %final off iter
    'FiniteDifferenceStepSize', 1e-9, ...
    'CheckGradients', true, ...
    'DiffMaxChange', 0.01, ...
    'OptimalityTolerance', 1e-9);

[vout, resnorm, residual, ~, ~, ~, jacobian] = lsqnonlin(ModelFunction, Start, Lower, Upper, options);


errors = get_errors(vout, residual, jacobian);

Coeff_array(N, :) = vout;
Errors_array(N, :) = errors;
end

%% look at single fit (last used)


eps1_model = Freal(vout, freq);
eps2_model = Fimag(vout, freq);

figure
hold on
plot(freq, eps2, '-xb')
plot(freq, eps2_model, '-r')
set(gca, 'xscale', 'log')

figure
hold on
plot(freq, eps1, '-xb')
plot(freq, eps1_model, '-r')
set(gca, 'xscale', 'log')

%% coeffs from current fit

figure
hold on
title('dEps')

errorbar(Temp, Coeff_array(:, 1), Errors_array(:, 1));
errorbar(Temp, Coeff_array(:, 4), Errors_array(:, 4));
errorbar(Temp, Coeff_array(:, 7), Errors_array(:, 7));
errorbar(Temp, Coeff_array(:, 10), Errors_array(:, 10));
set(gca, 'yscale', 'log')
xlim([255 355])
xlabel('T, K')
ylabel('dEps, 1')



figure
hold on
title('f_0')
errorbar(Temp, Coeff_array(:, 2), Errors_array(:, 2));
errorbar(Temp, Coeff_array(:, 5), Errors_array(:, 5));
errorbar(Temp, Coeff_array(:, 8), Errors_array(:, 8));
errorbar(Temp, Coeff_array(:, 11), Errors_array(:, 11));
set(gca, 'yscale', 'log')
xlim([255 355])
ylim([1e-6 1e12])
xlabel('T, K')
ylabel('f_0, Hz')



figure
hold on
title('alpha')
errorbar(Temp, Coeff_array(:, 3), Errors_array(:, 3));
errorbar(Temp, Coeff_array(:, 6), Errors_array(:, 6));
errorbar(Temp, Coeff_array(:, 9), Errors_array(:, 9));
errorbar(Temp, Coeff_array(:, 12), Errors_array(:, 12));
set(gca, 'yscale', 'log')
xlim([255 355])
ylim([0.01 1e1])
xlabel('T, K')
ylabel('alpha, 1')

yline(1)
yline(0)



%% plot all Eps1(f), Eps2(f)

figure('position', [510 144 609 804])

AX1 = subplot('position', [0.08    0.6698    0.875    0.292]);
hold on
set(gca, 'fontsize', 10)

AX1_res = subplot('position', [0.08     0.5330    0.875    0.102]);
set(gca, 'fontsize', 10)

AX2 = subplot('position', [0.08    0.2127    0.875    0.272]);
hold on
set(gca, 'fontsize', 10)

AX2_res = subplot('position', [0.08    0.0703    0.875    0.102]);
set(gca, 'fontsize', 10)

for N = 2:49

vout = Coeff_array(N, :);
% vout = coeff(N, :);

filename = ['./Data/freq/T_' num2str(N) '.TXT'];
[freq, eps1, eps2] = importfile(filename);

eps1_model = Freal(vout, freq);
eps2_model = Fimag(vout, freq);

XLIM = [5e-3 1e7];

axes(AX1)
cla
% plot(freq, eps1, 'xb')
errorbar(freq, eps1, eps1*relative_error, '.b')
plot(freq, eps1_model, '-r')
plot(freq, real(ColeCole(vout(1), vout(2), vout(3), freq)));
plot(freq, real(ColeCole(vout(4), vout(5), vout(6), freq)));
plot(freq, real(ColeCole(vout(7), vout(8), vout(9), freq)));
plot(freq, real(ColeCole(vout(10), vout(11), vout(12), freq)));
set(gca, 'xscale', 'log')
ylabel('eps1')
xlim(XLIM)
set(gca, 'XTick', [])
title(['T = ' num2str(Temp(N), '%3.0f') ' K'])

axes(AX1_res)
plot(freq, eps1-eps1_model, '-b')
set(gca, 'xscale', 'log')
xlim(XLIM)
ylim([-max(abs(ylim)) +max(abs(ylim))])
yline(0)

axes(AX2)
cla
% plot(freq, eps2, 'xb')
errorbar(freq, eps2, eps2*relative_error, '.b')
plot(freq, eps2_model, '-r')
plot(freq, -imag(ColeCole(vout(1), vout(2), vout(3), freq)));
plot(freq, -imag(ColeCole(vout(4), vout(5), vout(6), freq)));
plot(freq, -imag(ColeCole(vout(7), vout(8), vout(9), freq)));
plot(freq, -imag(ColeCole(vout(10), vout(11), vout(12), freq)));
set(gca, 'xscale', 'log')
ylabel('eps2')
xlim(XLIM)
set(gca, 'XTick', [])

axes(AX2_res)
plot(freq, eps2-eps2_model, '-b')
set(gca, 'xscale', 'log')
xlim(XLIM)
ylim([-max(abs(ylim)) +max(abs(ylim))])
yline(0)

drawnow
% pause(0.1)

% figure
% hold on
% plot(freq, eps1, '-xb')
% plot(freq, eps1_model, '-r')
% set(gca, 'xscale', 'log')



end




%% coeffs from fit.txt

figure
hold on
title('dEps')

plot(Temp, coeff(:, 1));
plot(Temp, coeff(:, 4));
plot(Temp, coeff(:, 7));
plot(Temp, coeff(:, 10));
set(gca, 'yscale', 'log')
xlim([255 355])
xlabel('T, K')
ylabel('dEps, 1')



figure
hold on
title('f_0')
plot(Temp, coeff(:, 2));
plot(Temp, coeff(:, 5));
plot(Temp, coeff(:, 8));
plot(Temp, coeff(:, 11));
set(gca, 'yscale', 'log')
xlim([255 355])
ylim([1e-6 1e12])
xlabel('T, K')
ylabel('f_0, Hz')



figure
hold on
title('alpha')
plot(Temp, coeff(:, 3));
plot(Temp, coeff(:, 6));
plot(Temp, coeff(:, 9));
plot(Temp, coeff(:, 12));
set(gca, 'yscale', 'log')
xlim([255 355])
ylim([0.01 1e1])
xlabel('T, K')
ylabel('alpha, 1')

yline(1)
yline(0)

