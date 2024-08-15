
vout(5) = 3;
Result = [];

N = 80;
Deps3_arr = linspace(vout(7)*0.2, vout(7)*5, N);
f3_arr = linspace(vout(8)*0.1, vout(8)*10, N);

for i = 1:numel(Deps3_arr)
    for j = 1:numel(f3_arr)
        Deps3 = Deps3_arr(i);
        f3 = f3_arr(j);

        Result(i, j) = sum(ModelFunction([vout(1:6) Deps3, f3, vout(9:12)]).^2);

    end
end

figure
imagesc(Deps3_arr, f3_arr, (Result'))
set(gca, 'Yscale', 'log')

xline(vout(7))
yline(vout(8))
caxis([1 40e3])




%%

Result = [];

N = 80
Deps2_arr = linspace(vout(4)*0.2, vout(4)*2, N);
f2_arr = linspace(vout(5)*0.05, vout(5)*2, N);

for i = 1:numel(Deps2_arr)
    for j = 1:numel(f2_arr)
        Deps2 = Deps2_arr(i);
        f2 = f2_arr(j);

        Result(i, j) = sum(ModelFunction([vout(1:3) Deps2, f2, vout(6:12)]).^2);

    end
end



figure
imagesc(Deps2_arr, f2_arr, log10(Result'))
set(gca, 'Yscale', 'log')

xline(vout(4))
yline(vout(5))
% caxis([1 40e3])

