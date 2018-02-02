function SR_LF = pb_lfsr_srcnn(LR_LF,mf,lf_name)

addpath('MATLAB/optical_flow/');
%--------------------------------------------------------------------------
% COMPUTE OPTICAL FLOW TO ALIGN THE LIGHT FIELD WITH CENTER VIEW
%--------------------------------------------------------------------------
if ~exist('FLOW/','dir')
    mkdir('FLOW/');
end

% Determine the filename where the flow vectors will be stored
forward_flow_filename = sprintf('FLOW/forward_%s.mat',lf_name);

if ~exist(forward_flow_filename,'file')
    % Compute the optical flow
    [u,v] = forward_optical_flow(LR_LF);
    % Save the optical flow vectors
    save(forward_flow_filename,'u','v');
else
    load(forward_flow_filename,'u','v');
end


% Use inverse warping to align all sub-aperture images to the centre view
LR_LF_align = inverse_warping(LR_LF,u,v);

clearvars LR_LF;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LIGHT FIELD DECOMPOSITION
%--------------------------------------------------------------------------

% Permute the channels of the aligned light field
LR_LF_align = permute(LR_LF_align,[1,2,4,5,3]);


%--- Display the progress of the sift-flow computation
msg = '  Matrix decomposition (R-Channel)...';
fprintf('%s',msg);
lengthLastMsg = length(msg);
pause(0.005);

[B_r,C_r] = matrix_decomposition(LR_LF_align(:,:,:,:,1));

%--- Clear the last entry
fprintf(repmat('\b', 1, lengthLastMsg));

%--- Display the progress of the sift-flow computation
msg = '  Matrix decomposition (G-Channel)...';
fprintf('%s',msg);
lengthLastMsg = length(msg);
pause(0.005);

[B_g,C_g] = matrix_decomposition(LR_LF_align(:,:,:,:,2));
%--- Clear the last entry
fprintf(repmat('\b', 1, lengthLastMsg));

%--- Display the progress of the sift-flow computation
msg = '  Matrix decomposition (B-Channel)...';
fprintf('%s',msg);
lengthLastMsg = length(msg);
pause(0.005);

[B_b,C_b] = matrix_decomposition(LR_LF_align(:,:,:,:,3));

%--- Clear the last entry
fprintf(repmat('\b', 1, lengthLastMsg));

% Encode the principal basis to get the principal bases as an image
[Ipb, param] = principal_basis_encoding(B_r(:,1), B_g(:,1), B_b(:,1),size(LR_LF_align));

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% SUPERRESOLUTION OF THE PRINCIPAL BASIS
%--------------------------------------------------------------------------

addpath('MATLAB\matconvnet\');
addpath('MATLAB\matconvnet\matlab\');
addpath('MATLAB\superresolution\srcnn\');

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
Ipb_sr = srcnn(uint8(Ipb*255),net);

% Convert the image in range [0,1]
Ipb_sr = double(Ipb_sr)/255;

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LIGHT FIELD RECONSTRUCTION
%--------------------------------------------------------------------------

% Decode the principal basis
[B_r(:,1), B_g(:,1),B_b(:,1)] = principal_basis_decoding(Ipb_sr,param);

% Initialize the SR_LF_align matrix
SR_LF_align = zeros(size(LR_LF_align));

% Reconstruct the super-ersolved aligned light field
SR_LF_align(:,:,:,:,1) = reshape(reshape(B_r*C_r,[size(LR_LF_align,1), ...
    size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
    [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
    size(LR_LF_align,4)]);
SR_LF_align(:,:,:,:,2) = reshape(reshape(B_g*C_g,[size(LR_LF_align,1), ...
    size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
    [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
    size(LR_LF_align,4)]);
SR_LF_align(:,:,:,:,3) = reshape(reshape(B_b*C_b,[size(LR_LF_align,1), ...
    size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
    [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
    size(LR_LF_align,4)]);

% % Permute the channels of the aligned light field
SR_LF_align = uint8(permute(SR_LF_align,[1,2,5,3,4]));

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% LIGHT FIELD RECONSTRUCTION
%--------------------------------------------------------------------------
% Restore the original disparities using forward warping
SR_LF = forward_warping(SR_LF_align,u,v);

% Specify the matrix completion method to be used here
matrix_completion_method = 'LMaFit';


% Inpaint the missing pixels using low-rank matrix completion
SR_LF = low_rank_matrix_completion(SR_LF,LR_LF, matrix_completion_method,lf_name);

function SR_LF = low_rank_matrix_completion(SR_LF, LR_LF,method,lf_name)

% Determine the filename where the flow vectors will be stored
inverse_flow_filename = sprintf('FLOW/inverse_%s.mat',lf_name);

% Compute sift flow to derive the flow vectors suitable to align the centre
% view to all the other views
if ~exist(inverse_flow_filename,'file')
    % Compute the optical flow
    [u,v] = inverse_optical_flow(LR_LF);
    % Save the optical flow vectors
    save(inverse_flow_filename,'u','v');
else
    load(inverse_flow_filename,'u','v');
end

% Extract the center views
It = SR_LF(:,:,:,5,5);

% Align the centre view to all sub-aperture images
k = 1; A = zeros(size(SR_LF,1),size(SR_LF,2),size(SR_LF,3),size(SR_LF,4)*size(SR_LF,5)-1);
for i = 1:size(SR_LF,4)
    for j = 1:size(SR_LF,5)
        if ~(i == 5 && j == 5)
            % Warp pixels from the sournce It to theIt aligned light field A
            for x = 1:size(It,2)
                for y = 1:size(It,1)
                    % Derive the new values of x and y
                    y_new = min(max(round(y - v(y,x,i,j)),1),size(It,1));
                    x_new = min(max(round(x - u(y,x,i,j)),1),size(It,2));
              
                    % Compute the warping
                    A(y,x,:,k) = It(y_new,x_new,:);
                end
            end
            k = k + 1;
        end
    end
end
%--------------------------------------------------------------------------
% Prepare the data matrix from the reconstructed light fields
M = reshape(SR_LF,[size(SR_LF,1)*size(SR_LF,2),size(SR_LF,3),size(SR_LF,4)*size(SR_LF,5)]);
M = permute(M,[1,3,2]);
% Prepare the auxiliary data which warps the centre view to all other views
A = reshape(A,[size(A,1)*size(A,2),size(A,3),size(A,4)]);
A = permute(A,[1,3,2]);
% Concatenate them to derive an augmented matrix
Maug = cat(2,M,A);
clearvars M A;
% Extract a matrix for each color channel
Raug = Maug(:,:,1); Gaug = Maug(:,:,2); Baug = Maug(:,:,3);
clearvars Maug;

% Derive the mask - holes are marked by 1
mask = Raug == -1;

% Compute the low-rank matrix completion to fill missing pixels
Raug = matrix_completion(Raug,method,mask,size(SR_LF));
Gaug = matrix_completion(Gaug,method,mask,size(SR_LF));
Baug = matrix_completion(Baug,method,mask,size(SR_LF));

% Reconstruct the Augmented matrix
Maug = cat(3,Raug,Gaug,Baug);
clearvars Raug Gaug Baug;

% Extrac the matrix of interest
M = Maug(:,1:81,:);
clearvars Maug;

% Permute the matrix
M = permute(M,[1,3,2]);

% Reconstruct the super-resolved light field
SR_LF = uint8(reshape(M,size(LR_LF)));



function X = matrix_completion(M,method,mask,sz)

figure(1);
imshow(uint8(reshape(M(:,1),[sz(1),sz(2)])));

% Derive the dimensions of the input matrix
[n1,n2] = size(M);

if strcmp(method,'LMaFit')
    addpath('MATLAB/matrix_completion/LMaFit/');
    addpath('MATLAB/matrix_completion/LMaFit/Utilities/');
    
    k = 10;

    %opts.rank_min = 10; % Minimum rank allowed
    Known = find(~mask);
    data = M(Known);
    opts.est_rank=2; %opts.est_rank=2;
    opts.rank_max=81;
    opts.Zfull=1;
    opts.tol=1e-5;
    opts.maxit=200;
    [B,C,~] = lmafit_mc_adp(n1,n2,k,Known,data,opts,sz);
    X = B*C;
end

% Replace all the known pixels in X
X(~mask) = M(~mask);





