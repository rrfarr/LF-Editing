function LF_super_resolution_analysis()
% This is the main function suitable to analyse the performance of the
% proposed light field super-resolution algorithm that computes existing
% super-resolution techniques on the principal basis. We decompose the
% light field into a set of basis elements X = BC and we super-resolve the
% principal basis B.

clc; close all;

addpath('MATLAB/light_field/');
addpath('MATLAB/superresolution/');

warning off;

%--------------------------------------------------------------------------
% SIMULATION CONFIGURATION
%--------------------------------------------------------------------------
sr_method   = 'bicubic';  % Super-resolution methods to consider. For now the 
                        % simulation will consider 'bicubic', 'srcnn' and 'vdsr'
                        
mf = 3;                 % Magnification factor

% Load the light fields to be processed in this simulation
[lf_names, datasets] = read_configuration('simulation.cfg');

% Determine the processing method
proc_method = get_sr_processing_method(sr_method);

if strcmp(sr_method,'srcnn') || strcmp(sr_method,'lf_srcnn')
    addpath('MATLAB/superresolution/srcnn/');
elseif strcmp(sr_method,'vdsr')
    addpath('MATLAB/superresolution/vdsr/');
end
% Determine the number of light fields to consider
N = size(lf_names,2);

fprintf('--------------------------------------------------------------\n');
if strcmp(proc_method,'independentSR') && (strcmp(sr_method,'srcnn') || (strcmp(sr_method,'lf_srcnn')))
    fprintf('Evaluating the performance of SRCNN\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/srcnn.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/srcnn/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/srcnn/',mf);
elseif strcmp(proc_method,'independentSR') && strcmp(sr_method,'bicubic')
    fprintf('Evaluating the performance of Bicubic interpolation\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/bicubic.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/bicubic/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/bicubic/',mf);
elseif strcmp(proc_method,'independentSR') && strcmp(sr_method,'vdsr')
    fprintf('Evaluating the performance of VDSR\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/vdsr.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/vdsr/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/vdsr/',mf);
end
fprintf('--------------------------------------------------------------\n');

if ~exist(out_img_foldername,'dir')
    mkdir(out_img_foldername);
end

if ~exist(out_LF_foldername,'dir')
    mkdir(out_LF_foldername);
end

% Open the file where the results will be stored
fid = fopen(out_filename,'w');

for n = 1:N
    % Get the dataset
    dataset = datasets{n};
    % Get the light field name
    lf_name = lf_names{n};
    
    out_img_filename = [out_img_foldername,lf_name,'.bmp'];
    out_LF_filename  = [out_LF_foldername,lf_name,'.mat'];
    
    % Derive the folder containing the light field
    dataset_foldername = sprintf('../LF-DATASET/%s/',dataset);
    
    % Load the high resolution light field
    if strcmp(dataset,'HCI')
        % Load the high resolution light field
        HR_LF = load_hci_lf(dataset_foldername, lf_name);
    elseif strcmp(dataset,'STANFORD')
        HR_LF = load_stanford_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'EPFL')
        HR_LF = load_epfl_lf(dataset_foldername,lf_name);
    elseif strcmp(dataset,'INRIA')
        HR_LF = load_inria_lf(dataset_foldername,lf_name);
    end
    
    % Permute the dimensions
    HR_LF = permute(HR_LF,[3,4,5,1,2]);
    
    % Generate the low-resolution light field
    LR_LF = lf_downsample(HR_LF,mf);
    
    if strcmp(proc_method,'independentSR')
        % The light field will be restored one sub-aperture image at a time
        % and will not exploit the light field structure
        if strcmp(sr_method,'bicubic')
            SR_LF = lf_bicubic(LR_LF);
        elseif strcmp(sr_method,'srcnn') || strcmp(sr_method,'lf_srcnn')
            SR_LF = lf_srcnn(LR_LF,mf);
        elseif strcmp(sr_method,'vdsr')
            SR_LF = lf_vdsr(LR_LF);
        end
            
    elseif strcmp(proc_method,'coherentSR')
    end
    
    % Extract the centre view
    Ic = SR_LF(:,:,:,5,5);
    
    % Write the center view
    imwrite(Ic,out_img_filename);
    
    % Compute the psnr evaluation
    psnr_val = quality_analysis(HR_LF, SR_LF);
    
    fprintf('%s, %s, %0.4f\n',lf_name,dataset,psnr_val);
    fprintf(fid,'%s, %s, %0.4f\n',lf_name,dataset,psnr_val);
    
    % Save the restored light field
    save(out_LF_filename,'SR_LF');
end
fclose('all');

