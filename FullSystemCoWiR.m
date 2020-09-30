%% COWIR Radar System
%{

    Sean Holloway
    CoWiR (Compact Wind Radar) System
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.

    README:
    Place all .bin files to be parsed in the 'Input Data' folder.
    Run as normal to process all files in folder.

    Uncomment lines 46 and 47 to replace with specific file names.

    Change line 52 from "%" to "%{" to skip parsing if previously
    completed.


%}


%% Housekeeping
clear variables
close all
addpath(genpath('MAT Files'));
tic

%% User Options

subtract_signal = true;
% subtract_signal = false;

%% Setup

% Discover all files in 'Input Data' folder
fold = dir('Input Data');
k = 1;
for i = 1:length(fold)
    if not(fold(i).isdir)
        files{k} = fold(i).name(1:end-4);
        k = k+1;
    end
end

% Overwrite file names
files = {'beach_straight_out_1'};


%% Interpret Raw Data from Files

%
for file_loop = 1:length(files)
    
    seq_name = files{file_loop};
    if ~isfile(['Parsed Data\', seq_name, '_1.mat'])
        message = ['Parsing data for file ', seq_name, sprintf('"\nFile # %d of %d.', file_loop, length(files))];
        disp(message);
        toc
        ParseFiles
    else
        message = ['File "', seq_name, '" has already been parsed.'];
        disp(message);
    end
end
%}

%% Signal & Data Processing

for file_loop = 1:length(files)
    
    seq_name = files{file_loop};    
    message = ['Processing radar cube for file ', seq_name, sprintf('"\nFile # %d of %d.', file_loop, length(files))];
    disp(message);
    toc
    
    for chan_num = 1:3
        
        message = ['Processing channel ', sprintf('%d',chan_num)];
        disp(message);
    
        numfiles = length(dir(['Parsed Data\', seq_name, '_*.mat']));
        
        for seqnumber = 1:numfiles
            
            % Save settings to pass into signal processing
            save('MAT Files/setup.mat', 'seq_name', 'seqnumber')
            
            % Perform R-D FFTs and binary integration
            RadarSignalProcessing_v1
            
            % Save each radar cube into full sequence cube
            seq_cube(:,:,seqnumber) = integrated_signal;
            save('MAT Files/seq_cube.mat', 'seq_cube', 'range_axis', 'doppler_axis', '-v7.3');
            
        end
        
        % Process R-D data
        RadarDataProcessing_v1
    end
end

%% Wind Velocity Estimation

for file_loop = 1:length(files)
    
    seq_name = files{file_loop};
    
    message = ['Processing velocity data for file ', seq_name, sprintf('"\nFile # %d of %d.', file_loop, length(files))];
    disp(message);
    toc
    
    % Estimate True Wind Velocity
    WindVelocity
end


