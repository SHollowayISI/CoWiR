%% Radar Signal Processing for CoWiR Project
%{
    
    Sean Holloway
    2/10/2020
    Version 1
    Radar signal processing for CoWiR project.

    Range-Doppler post-processing for recorded signal from ParseFiles.m
    results.

    Working as of 2/10/2020

%}

%% Housekeeping
% clear variables
clear xlabel;
clear ylabel;
close all;
% tic

addpath(genpath('Parsed Data'));
addpath(genpath('MAT Files'));
addpath(genpath('Functions'));

c = physconst('LightSpeed');

%% Variables

% Assigned Variables
fc = 14.4e9;                    % Operating frequency in Hz
fadc = 100e6;                   % Sample frequency in Hz
fd = 20e6;                      % Sample frequency after decimation
tm = 10e-6;                     % Sweep time in seconds
bw = 50e6;                      % Sweep bandwidth in Hz

adc_samples = 1250;             % Number of samples recorded in a chirp
num_samples = 250;              % Number of samples after decimation
samples_per_chirp = 212;        % Number of samples to keep
range_fft_size = 256;           % Number of samples to use in range FFT

% NOTE: Currently using numbers from lab setup
adc_samples = 1000;%1250;             % Number of samples recorded in a chirp
num_samples = 200;%250              % Number of samples after decimation
samples_per_chirp = 170;%220;        % Number of samples to keep
range_fft_size = 256;           % Number of samples to use in range FFT

num_chirps = 1024;              % Number of chirps per frame
num_frames = 10;                % Number of frames per CPI
num_cpi = 1;                   % Number of recorded CPI


% Derived Variables
lambda = c/fc;
ts = 1/fd;
dec_factor = adc_samples/num_samples;
drop_samples = num_samples - samples_per_chirp;
sweep_slope = bw/tm;
chirp_rate = 1/tm;
frame_time = tm*num_chirps;

range_res = (num_samples/samples_per_chirp)*c/(2*bw);
vel_res = lambda/(2*frame_time);

range_res = 2.34;%3.02;%
vel_res = 0.156;


%% Load signal from file

% Set up file name
if isfile('MAT Files/setup.mat')
    load('setup.mat');
else
%     seq_name = 'out_to_sea4';
%     seqnumber = 1;
end

filename = sprintf([seq_name, '_%d.mat'], seqnumber);

% Load data
data_in = load(filename);
loadchan = sprintf('signal = data_in.chan%d;', chan_num);
eval(loadchan);
sync = data_in.sync;

% sync = cowirwalkingtoward(1:10240000,3);
% signal = cowirwalkingtoward(1:10240000,8) + 1i*cowirwalkingtoward(1:10240000,6);

%% Signal shaping and decimation

% Reshape and isolate IQ signal
% signal = data_in(:,2) + 1i*data_in(:,1);
% sync = data_in(:,4);

% Remove samples before first sync pulse and after final sync pulse
[~, loc] = findpeaks(sync, 'MinPeakDistance', 10, 'MinPeakProminence', 100);
% loc_dist = (loc(end)-loc(1));
% signal = padarray(signal(loc(1):(loc(end)-1)), (10240000 - loc_dist), 0, 'post');

signal = signal(loc(1):(loc(end)-1));

% Decimate signal by factor of 5
signal = decimate(signal, dec_factor);

% Reshape to fast-time x slow-time
signal = reshape(signal, num_samples, []);

% Remove mean signal from total signal
if subtract_signal
    signal = signal - mean(signal, 2);
end

% Remove early samples from each chirp and zero-pad to correct lengths
signal = signal((drop_samples+1):end, :);

range_pad = range_fft_size - samples_per_chirp;
doppler_pad = num_chirps*num_frames*num_cpi - size(signal,2);

signal = padarray(signal, [range_pad doppler_pad], 0, 'post');

% Reshape to fast-time x slow-time x frame x CPI
signal = reshape(signal, range_fft_size, num_chirps, num_frames, num_cpi);

%% Range FFT & Windowing

% Conjugate signal to convert negative frequency to positive
signal = conj(signal);

% Apply Hanning window to fast time dimension
h = hann(range_fft_size);
signal = h.*signal;

% Create range axis
k_r = (1:(range_fft_size/2)) - 0.5;
range_axis = k_r * range_res;

% Apply range FFT
signal = fft(signal,range_fft_size,1);

% Drop second half of range FFT frequencies
signal = signal(1:(end/2), :, :, :);

%% Doppler FFT & Windowing

% Apply Blackman-Harris windowing
h = blackmanharris(num_chirps);
signal = h'.*signal;

% Create axis
k_d = (-num_chirps/2:num_chirps/2);
doppler_axis = k_d * vel_res;

% Perform FFT and shift
signal = fft(signal, num_chirps, 2);
signal = fftshift(signal, 2);

% Wrap end of FFT shfit
signal(:,(end+1),:,:) = signal(:,1,:,:);

%% Non-Coherent Integration

% Calculate power of signal
power_signal = abs(signal).^2;

% Calculate mean power by frame
integrated_signal = squeeze(mean(power_signal, 3));

%% Plot Range-Doppler Intensity Heatmap
%
figure('Name', 'Range-Doppler Intensity Heatmap')
surfc(doppler_axis, range_axis, 10*log10(integrated_signal), 'EdgeColor', 'none')
set(gca,'YDir','normal')
xlabel('Velocity [m/s]')
ylabel('Range [m]')
xlim([-20 20])
ylim([0 100])
zlim([60 120])
%}

% toc;

%% Save full cube to file
%
% disp('Processing complete, saving file...')
filename = 'MAT Files/processedCube.mat';
save(filename,'integrated_signal', 'range_axis', 'doppler_axis', '-v7.3');
% toc;
% disp(['Signal saved as "', filename, '"'])
%}



