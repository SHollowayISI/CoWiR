%% File Parsing for CoWiR Project
%{
    
    Sean Holloway
    2/10/2020
    Version 1
    File parsing for CoWiR project.

    Converts data from .txt files to usable .mat matrix file.

    Working as of 2/10/2020

%}

%% Housekeeping

% clear variables;
close all;
% tic

addpath(genpath('Parsed Data'));
addpath(genpath('Input Data'));

%% Load .txt files and convert to .mat arrays
%{
% seqnumber = 1;
% seqname = 'd';

for seqnumber = 1:10
    filename = sprintf(['cowir_walking_towardseq_matrix_s%d.txt'] ,seqnumber);
    fileID = fopen(filename, 'r');
    data_in = reshape(fscanf(fileID, '%f,%f,%f,%f'), 4, [])';
    
    filename = sprintf(['Parsed Data/cowir_walking_toward%d.mat'], seqnumber);
    save(filename, 'data_in');
    toc;
    disp(['Signals saved as "', filename, '"'])
    
end
%}

%% Load .bin files and convert to to .mat arrays

% seq_num = 4;
% seq_name = 'out_to_sea';

read_filename = [seq_name, '.bin'];

fileID=fopen(read_filename,'r');
outbuff=fread(fileID,'int16');

size_buff = length(outbuff);
outbuff = padarray(outbuff, 8-mod(size_buff,8), 0, 'post');

outbuff = reshape(outbuff, 8, [])';
[size_buff,dimx]=size(outbuff);

numfiles = ceil(size_buff/10240000);
outbuff = padarray(outbuff, numfiles*10240000-size(outbuff,1), 0, 'post');
for i = 1:numfiles
    
    ind_range = ((i-1)*10240000 + 1):(i*10240000);
    
    chan1 = outbuff(ind_range, 7) + 1i*outbuff(ind_range, 5);
    chan2 = outbuff(ind_range, 4) + 1i*outbuff(ind_range, 2);
    chan3 = outbuff(ind_range, 8) + 1i*outbuff(ind_range, 1);
    sync = outbuff(ind_range, 3);
    
    write_filename = ['Parsed Data/', seq_name, sprintf('_%d.mat', i)];
    save(write_filename, 'chan1', 'chan2', 'chan3', 'sync');,
    
end






