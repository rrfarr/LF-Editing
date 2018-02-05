function quality_performance_analysis_SR()

clc; close all; clear all;

%--------------------------------------------------------------------------
% Configure the analysis
%--------------------------------------------------------------------------
addpath('../MATLAB/general/');

mf = 3; % magnification factor

%--------------------------------------------------------------------------

% Specify the foldername where the results are stored
results_foldername = sprintf('../RESULTS/superresolution/x%d/',mf);

% Get the file list
list = get_file_list(results_foldername,'csv');

% Determine the number of files to be considered
Nfiles = size(list,1);

header{1} = 'LF-name';
header{2} = 'Dataset';
for n = 1:Nfiles
    % Derive the method considered
    method = list{n}(1:end-4);
    
    if strcmp(method,'bicubic')
        method  = 'Bicubic';
        cell_id = 1;
    elseif strcmp(method,'bm_pca_rr')
        method  = 'BM+PCA+RR';
        cell_id = 3;
    elseif strcmp(method,'pca_rr')
        method  = 'PCA+RR';
        cell_id = 2;
    else
        error('This method is not considered\n');
    end
    
    % Open the file
    fid = fopen([results_foldername, list{n}],'r');
    
    % Put the method as a header of the xls file
    header{2+cell_id} = method;
    % Read the content of this file
    [lf_name, ds_name, psnr, ssim] = csv_read_SR_quality_file(fid);
    
    % Close the file
    fclose(fid);
end

% Derive the xls filename
xls_foldername = sprintf('../RESULTS/superresolution/x%d/performance_analysis/',mf);

if ~exist(xls_foldername,'dir')
    mkdir(xls_foldername);
end

% Creat the xls filename
xls_filename = [xls_foldername,'quality_analysis.xls'];


function [lf_name, ds_name, psnr, ssim] = csv_read_SR_quality_file(fid)

k = 1;
while ~feof(fid)
    % Read a line from the csv file
    line = fgetl(fid);
    % Derive the index of comma
    idx_comma = strfind(line,',');
    lf_name{k} = line(1:idx_comma(1)-1);
    ds_name{k} = line(idx_comma(1)+1:idx_comma(2)-1);
    psnr(k)    = str2double(line(idx_comma(2)+1:idx_comma(3)-1));
    ssim(k)    = str2double(line(idx_comma(3)+1:end));
    k = k + 1;
end




