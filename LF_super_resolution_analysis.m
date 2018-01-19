function LF_super_resolution_analysis()
% This is the main function suitable to analyse the performance of the
% proposed light field super-resolution algorithm that computes existing
% super-resolution techniques on the principal basis. We decompose the
% light field into a set of basis elements X = BC and we super-resolve the
% principal basis B.

clc; close all;

addpath('MATLAB/light_field/');
warning off;

%--------------------------------------------------------------------------
% SIMULATION CONFIGURATION
%--------------------------------------------------------------------------
sr_method   = 'vdsr';  % Super-resolution methods to consider. For now the 
                        % simulation will consider 'bicubic', 'srcnn' and 'vdsr'
                        
proc_method = 1;        % This specifies if we are going to compute super
% proc_method = 2       % resolution on each image separately (1) or else 
                        % restoring the principal basis only (2)
                        
mf = 3;                 % Magnification factor

lf_names = {'Bikes', ...
    'Bench_in_Paris' ...
    'Fountain_&_Bench', ...
    'Friends_1', ...
    'Sphynx', ...
    'Bee_2', ...
    'Duck', ...
    'Fruits', ...
    'Rose', ...
    'Mini', ...
    'Chess',...
    'Bunny',...
    'Lego Bulldozer', ...
    'Lego Truck', ...
    'Lego Knights'
    };
datasets = {'EPFL', ...
    'EPFL', ...
    'EPFL', ...
    'EPFL', ...
    'EPFL', ...
    'INRIA', ...
    'INRIA', ...
    'INRIA', ...
    'INRIA', ...
    'INRIA', ...
    'STANFORD', ...
    'STANFORD', ...
    'STANFORD', ...
    'STANFORD', ...
    'STANFORD'
    };
if strcmp(sr_method,'srcnn')
    addpath('MATLAB/superresolution/srcnn/');
elseif strcmp(sr_method,'vdsr')
    addpath('MATLAB/superresolution/vdsr/');
end
% Determine the number of light fields to consider
N = size(lf_names,2);

fprintf('--------------------------------------------------------------\n');
if proc_method == 1 && strcmp(sr_method,'srcnn')
    fprintf('Evaluating the performance of SRCNN\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/srcnn.csv',mf);
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/srcnn/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/srcnn/',mf);
elseif proc_method == 1 && strcmp(sr_method,'bicubic')
    fprintf('Evaluating the performance of Bicubic interpolation\n');
    out_filename = sprintf('RESULTS/superresolution/x%d/bicubic.csv',mf);;
    out_img_foldername = sprintf('RESULTS/superresolution/x%d/centre_view/bicubic/',mf);
    out_LF_foldername  = sprintf('RESULTS/superresolution/x%d/LF/bicubic/',mf);
elseif proc_method == 1 && strcmp(sr_method,'vdsr')
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
    
    % Initialize the super-resolved light field
    SR_LF = zeros(size(LR_LF),'uint8');
    if proc_method == 1
        % The light field will be restored one sub-aperture image at a time
        % and will not exploit the light field structure
        for u = 1:size(LR_LF,4)
            for v = 1:size(LR_LF,5)
                % Extract the current sub-aperture image
                Ilr = LR_LF(:,:,:,u,v);
                if strcmp(sr_method,'bicubic')
                    Isr = Ilr;
                elseif strcmp(sr_method,'srcnn')
                    addpath('MATLAB\matconvnet\');
                    addpath('MATLAB\matconvnet\matlab\');

                    %run matconvnet/matlab/vl_setupnn;
                    run vl_setupnn

                    model_filename = sprintf('MATLAB\\superresolution\\srcnn\\mf_%d_model.mat',mf);
                    load(model_filename);

                    % Test that the convolution function works properly
                    try
                        vl_nnconv(single(1),single(1),[]) ;
                    catch
                        warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
                        vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
                            'enableImreadJpeg', false) ;
                    end

                    % Super-resolve the current sub-aperture image
                    Isr = srcnn(Ilr,net);
                elseif strcmp(sr_method,'vdsr')
                    addpath('MATLAB\matconvnet\');
                    addpath('MATLAB\matconvnet\matlab\');

                    %run matconvnet/matlab/vl_setupnn;
                    run vl_setupnn

                    addpath('MATLAB\superresolution\vdsr\utils\');
                    load('MATLAB\superresolution\vdsr\VDSR_Official.mat');
                    
                    % Test that the convolution function works properly
                    try
                        vl_nnconv(single(1),single(1),[]) ;
                    catch
                        warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
                        vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
                            'enableImreadJpeg', false) ;
                    end
                    
                    % Super-resolve the current sub-aperture image
                    Isr = vdsr(Ilr,model,mf);

                end
                % Put the restored sub-aperture image in the restored light
                % field
                SR_LF(:,:,:,u,v) = Isr;
            end
        end
    elseif proc_method == 2
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


