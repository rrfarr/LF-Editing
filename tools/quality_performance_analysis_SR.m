function quality_performance_analysis_SR()

clc; close all; clear all;
%--------------------------------------------------------------------------
% Configure the analysis
%--------------------------------------------------------------------------
addpath('MATLAB/');
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
    elseif strcmp(method,'pblfsr_srcnn')
        method = 'Proposed';
        cell_id = 6;
    elseif strcmp(method,'srcnn')
        method = 'LF-SRCNN';
        cell_id = 4;
    elseif strcmp(method,'vdsr')
        method = 'LF-VDSR';
        cell_id = 5;
    else
        error('This method is not considered\n');
    end
    
    % Open the file
    fid = fopen([results_foldername, list{n}],'r');
    
    % Put the method as a header of the xls file
    header{2+cell_id} = method;
    % Read the content of this file
    [lf_name, ds_name, psnr, ssim] = csv_read_SR_quality_file(fid);

    psnr_cell(:,cell_id) = psnr;
    ssim_cell(:,cell_id) = ssim;
    
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

xlswrite(xls_filename,header,'PSNR');
xlswrite(xls_filename,lf_name,'PSNR',sprintf('A%d:A%d',2,1+length(lf_name)));
xlswrite(xls_filename,ds_name,'PSNR',sprintf('B%d:B%d',2,1+length(lf_name)));
xlswrite(xls_filename,psnr_cell,'PSNR',sprintf('C%d:%s%d',2,'C'+size(psnr_cell,2)-1,1+length(lf_name)));

xlswrite(xls_filename,header,'SSIM');
xlswrite(xls_filename,lf_name,'SSIM',sprintf('A%d:A%d',2,1+length(lf_name)));
xlswrite(xls_filename,ds_name,'SSIM',sprintf('B%d:B%d',2,1+length(lf_name)));
xlswrite(xls_filename,ssim_cell,'SSIM',sprintf('C%d:%s%d',2,'C'+size(psnr_cell,2)-1,1+length(lf_name)));

