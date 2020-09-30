%% Radar Data Processing for CoWiR Project
%{
    
    Sean Holloway
    2/10/2020
    Version 1
    Radar data processing for CoWiR project.

    Data processing using results of RadarSignalProcessing_v1.

    Working as of 2/10/2020

%}

%% Housekeeping
% clear variables;
close all;

addpath(genpath('MAT Files'));


%% Load signal from file

% Load data
% if isfile('MAT Files/seq_cube.mat')
%     cube_in = load('MAT Files/seq_cube.mat');
%     
%     % Unload .mat file
%     range_axis = cube_in.range_axis;
%     doppler_axis = cube_in.doppler_axis;
%     cube_in = 10*log10(cube_in.seq_cube);
%     
% else
    filename = 'seq_cube.mat';
    cube_in = load(filename);
    
    % Unload .mat file
    range_axis = cube_in.range_axis;
    doppler_axis = cube_in.doppler_axis;
    
    % Log import:
    cube_in = 10*log10(cube_in.seq_cube);
    % Remove infinities
    cube_in(cube_in == -Inf) = -10000;

%     % Linear import:
% 	cube_in = cube_in.seq_cube;
    
% end


%% Velocity determination

% Set range-velocity area of focus
max_range = 60;
max_vel = 40;
[~, max_bin] = min(abs(range_axis - max_range));
max_bin = max_bin - 1;

velocity_est = zeros(max_bin, size(cube_in,3));

% Choose CPI to view
for cpi_num = 1:size(cube_in,3)
    
    cube = cube_in(:,:,cpi_num);
    
    fit_cube = cube - median(cube, 2);
    
    %% Leakage signal removal via interpolation
    
    % Linear interpolation between center slices
    %
    omitwidth = 5;
    int_cube(:,:,cpi_num) = cube;
    int_cube(:,ceil(end/2)-omitwidth:ceil(end/2)+omitwidth,cpi_num) = ...
        ones(size(cube,1),(2*omitwidth+1)).*linspace(0,1,(omitwidth*2 + 1)).*...
        (cube(:,ceil(end/2)+omitwidth) - cube(:,ceil(end/2)-omitwidth))...
        + cube(:,ceil(end/2)-omitwidth);
    %}
    


    
    
    %% Velocity estimation via curve fitting
    
    % Bi-Gaussian fitting for log plot across Doppler bins
    %{
for i = 1:max_bin
    
    f = fit(doppler_axis.', fit_cube(i,:)', 'gauss2', ...
        'Lower', [0, -1, 0, 0, -max_vel, 0], ...
        'Upper', [Inf, 1, max_vel/20, Inf, max_vel, max_vel/4]);
    
    fit_cube(i,:) = fit_cube(i,:) - f.a1*exp(-((doppler_axis-f.b1)/f.c1).^2);
    
    fit_cube(i,:) = movmean(fit_cube(i,:),10);
    
    g = fit(doppler_axis', fit_cube(i,:)', 'gauss1', ...
        'Lower', [0, -max_vel, 0], ...
        'Upper', [Inf, max_vel, max_vel/4]);
    
    velocity_est(i) = g.b1;
end
    %}
    
    % Single Gaussian fitting of linear interpolated data
    %
    for i = 1:max_bin
        
        f{i}{cpi_num} = fit(doppler_axis.', (int_cube(i,:,cpi_num)-median(int_cube(i,:,cpi_num)))', ...
            'gauss1', ...
            'Lower', [0, -max_vel, 0], ...
            'Upper', [Inf, max_vel, max_vel/4]);
        
        velocity_est(i,cpi_num) = f{i}{cpi_num}.b1;
        
    end
    %}
    
end


%% Visualization

view_cpi = 1;

% Plot range slices
%
figure('Name', 'Range Slices')
for i = 1:20
%     [~, idx] = min(abs(range_axis - 5*i));
    idx = i;
    range_bin = range_axis(idx);
    plot(doppler_axis, int_cube(idx,:,view_cpi), ...
        'DisplayName', sprintf('%0.1f meters', range_axis(i)), ...
        'LineWidth', 1.5);
    set(gca, 'FontWeight', 'bold', 'XLim', [-20, 20]);
    xlabel('Velocity [m/s]');
    ylabel('FFT Log Intensity');
    hold on
end
title('Velocity Response at Different Ranges')
legend
%}

% Plot range slices without legend, for export reasons
%{
figure('Name', 'Range Slices No Legend')
for i = 1:20
    [~, idx] = min(abs(range_axis - 5*i));
    range_bin = range_axis(idx);
    plot(doppler_axis, int_cube(:,idx), ...
        'DisplayName', sprintf('%d meters', 5*i));
    set(gca, 'XLim', [-10, 10]);
    xlabel('Velocity [m/s]');
    ylabel('FFT Log Intensity');
    hold on
end
%}

% Plot range slices with Gaussian fit
%{
figure('Name', 'Gaussian Fit')
range_bin = 5;
plot(doppler_axis, int_cube(range_bin,:,view_cpi)-median(int_cube(range_bin,:,view_cpi)),...
    'LineWidth', 1, ...
    'DisplayName', sprintf('Range Slice at %0.0fm', 5*round(range_axis(range_bin)/5)))
legend
hold on
plot(f{range_bin}{view_cpi})
xlim([-max_vel, max_vel])
xlabel('Velocity [m/s]')
ylabel('Adjusted Log Intensity')
set(gca, 'FontWeight', 'bold');
title('Gaussian Fit of Velocity Response')
%}

% Plot Gaussian-fit velocity estimates
%
figure('Name', 'Velocity Estimate')
plot(range_axis(1:max_bin), velocity_est)
grid on;
set(gca, 'YLim', [-10, 10], 'FontWeight', 'bold');
xlabel('Range [m]');
ylabel('Estimated Velocity [m/s]');
hold on

mean_vel = mean(velocity_est,2);
std_vel = std(velocity_est')';
plot(range_axis(1:max_bin), mean_vel, ...
    'LineWidth', 2, 'Color', 'k')

% Add legend
for ind = 1:size(velocity_est,2)
    data_titles{ind} = sprintf('Frame #%d', ind);
end
data_titles{ind+1} = 'Mean';
legend(data_titles)
title('Velocity Estimate per Frame')

% Plot mean and stdev of above plot
figure('Name', 'Mean Velocity and Standard Deviation')
plot(range_axis(1:max_bin), mean_vel, ...
    'LineWidth', 2, 'Color', 'k', ...
    'DisplayName', 'Mean')
hold on
plot(range_axis(1:max_bin), mean_vel-std_vel, ...
    range_axis(1:max_bin), mean_vel+std_vel, ...
    'LineWidth', 1, 'Color', 'r', ...
    'DisplayName', '+/-1 Dev')
hold on
plot(range_axis(1:max_bin), mean_vel-2*std_vel, ...
    range_axis(1:max_bin), mean_vel+2*std_vel, ...
    'LineWidth', 0.5, 'Color', '#EDB120', ...
    'DisplayName', '+/-2 Dev')
grid on;
legend;
set(gca, 'YLim', [-10, 10], 'FontWeight', 'bold');
xlabel('Range [m]');
ylabel('Estimated Velocity [m/s]');
%}

% Plot Range-Doppler surface
%
figure('Name', 'Range-Doppler Intensity Heatmap')
surfc(doppler_axis, range_axis, int_cube(:,:,view_cpi), 'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Velocity [m/s]')
ylabel('Range [m]')
xlim([-20 20])
ylim([0 100])
zlim([60 120])
%}

%% Save plots to file
%
filepath = ['Figures\', ...
    seq_name, '\'];
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

FolderName = filepath;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    savefig(FigHandle, fullfile(FolderName, [seq_name, '_ch', sprintf('%d', chan_num), '_', FigName, '.fig']));
    saveas(FigHandle, fullfile(FolderName, [seq_name, '_ch', sprintf('%d', chan_num), '_', FigName, '.png']));
end

% close all
%}

%% Save velocity data

filepath = 'MAT Files\';
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

filename = [filepath, seq_name, '_ch', sprintf('%d', chan_num), '.mat'];

save(filename, 'range_axis', 'mean_vel');





