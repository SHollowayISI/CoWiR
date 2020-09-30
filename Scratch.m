close all;

test_slice = squeeze(int_cube(3,:,1));

test_slice = test_slice-median(test_slice);

% test_slice_mod = movmean(test_slice,10);
test_slice_mod = test_slice;

f = fit(doppler_axis', test_slice_mod', 'gauss2', ...
    'Lower', [0, -1, 0, 0, -max_vel, 0], ...
    'Upper', [Inf, 1, max_vel/20, Inf, max_vel, max_vel/4]);

test_slice_sub = test_slice_mod - f.a1*exp(-((doppler_axis-f.b1)/f.c1).^2);

test_slice_sub = movmean(test_slice_sub,10);

g = fit(doppler_axis', test_slice_sub', 'gauss1', ...
    'Lower', [0, -max_vel, 0], ...
    'Upper', [Inf, max_vel, max_vel/4]);

figure;
plot(f)
hold on
plot(g)
hold on
plot(doppler_axis, test_slice)
hold on
plot(doppler_axis, test_slice_mod)
hold on
plot(doppler_axis, test_slice_sub)
