function SR_LF = vdsr_LFSR(LR_LF,sz,mf)
% This function is a wrapper script that allows us to use the code provided
% by the authors of the paper Chao Dong, Chen Change Loy, Kaiming He, Xiaoou 
% Tang. Image Super-Resolution Using Deep Convolutional Networks, IEEE 
% Transactions on Pattern Analysis and Machine Intelligence (TPAMI), 2015.
% The SRCNN was also used in a recent paper in lightfield super-resolution
% as mentioned by the same authors in Y. Yoon, H. G. Jeon, D. Yoo, J. Y. Lee
% and I. S. Kweon, "Learning a Deep Convolutional Network for Light-Field 
% Image Super-Resolution," 2015 IEEE International Conference on Computer 
% Vision Workshop (ICCVW), Santiago, 2015, pp. 57-65.
% This method receives the low-resolution lightfield LR_LF, and computes
% super-resolution so that each sub-aparture images has a dimension sz with
% a magnification factor provided in scale.
%
% Input: LR_LF - low-resolution lightfield
%        sz    - dimensions of each sub-aparture image
%        scale - magnification factor
%
% Output: SR_LF - the super-resolved lightfield image.
%
% Reuben Farrugia
% 30/07/2016

addpath('MATLAB\matconvnet\');
addpath('MATLAB\matconvnet\matlab\');

%run matconvnet/matlab/vl_setupnn;
run vl_setupnn

addpath('MATLAB\vdsr\utils\');
load('MATLAB\vdsr\VDSR_Official.mat');

% Test that the convolution function works properly
try
  vl_nnconv(single(1),single(1),[]) ;
catch
  warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
  vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
               'enableImreadJpeg', false) ;
end


%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
%load the pre-computed model (one for each scale)
%fprintf('------------------------------------------------------------------\n');
%fprintf('Load the pre-computed SRCNN model\n');

% Initialize super-resolved lightfield
SR_LF = zeros(sz(1),sz(2),size(LR_LF,3));

% Determine the number of lightfields
N = size(SR_LF,3);

%--------------------------------------------------------------------------
% Compute the lightfield super-resolution using SRCNN method
%--------------------------------------------------------------------------
% Take the time when the simulation started
for n = 1:N
    % Load the low-resolution image
    I_LR = uint8(LR_LF(:,:,n));
    
    % Initialize the timer
    img_strt = tic;
    
    % Upscale the low-resolution image
    I_LR = imresize(I_LR,sz,'bicubic');
    
    % Convert the low-resolution image to double
    I_LR = double(I_LR)/255;
    
    % Compute the super-resolution
    im = VDSR_Matconvnet(I_LR,model,mf,0);
    im = im * 255;
    %---------------------------------------------------------------------
    
    % Number of seconds to process sub aparture image
    img_exp = toc(img_strt);
    
    % Derive the total time 
    tot_time = N*img_exp;
    
    % Derive the total time already processed
    tot_proc = n*img_exp;
    
    % Derive the total time left
    tot_left = tot_time - tot_proc;
    
    % Derive the number of seconds left
%    left = seconds2human(tot_left);
    
    % Derive the super-resolved image
    I_SR   = uint8(im);
    
    % Convert the image to uint8
    SR_LF(:,:,n) = I_SR;
    
    % Report the number of sub-aparture images processed
%    fprintf('%d out of %d sub-aparture imgs ready (%s left)\n',n,N,left);
%    fflush(stdout);
end
SR_LF = uint8(SR_LF);
%fprintf('------------------------------------------------------------------\n');
%fflush(stdout);

