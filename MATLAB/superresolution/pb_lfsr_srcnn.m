function SR_LF = pb_lfsr_srcnn(LR_LF,mf,lf_name)

% addpath('MATLAB/optical_flow/');
% %--------------------------------------------------------------------------
% % COMPUTE OPTICAL FLOW TO ALIGN THE LIGHT FIELD WITH CENTER VIEW
% %--------------------------------------------------------------------------
% if ~exist('FLOW/','dir')
%     mkdir('FLOW/');
% end
% 
% % Determine the filename where the flow vectors will be stored
% flow_filename = sprintf('FLOW/%s.mat',lf_name);
% 
% if ~exist(flow_filename,'file')
%     % Compute the optical flow
%     [u,v] = optical_flow(LR_LF);
%     % Save the optical flow vectors
%     save(flow_filename,'u','v');
% else
%     load(flow_filename,'u','v');
% end
% 
% % Use inverse warping to align all sub-aperture images to the centre view
% LR_LF_align = inverse_warping(LR_LF,u,v);
% 
% clearvars LR_LF;
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % LIGHT FIELD DECOMPOSITION
% %--------------------------------------------------------------------------
% 
% % Permute the channels of the aligned light field
% LR_LF_align = permute(LR_LF_align,[1,2,4,5,3]);
% 
% 
% %--- Display the progress of the sift-flow computation
% msg = '  Matrix decomposition (R-Channel)...';
% fprintf('%s',msg);
% lengthLastMsg = length(msg);
% pause(0.005);
% 
% [B_r,C_r] = matrix_decomposition(LR_LF_align(:,:,:,:,1));
% 
% %--- Clear the last entry
% fprintf(repmat('\b', 1, lengthLastMsg));
% 
% %--- Display the progress of the sift-flow computation
% msg = '  Matrix decomposition (G-Channel)...';
% fprintf('%s',msg);
% lengthLastMsg = length(msg);
% pause(0.005);
% 
% [B_g,C_g] = matrix_decomposition(LR_LF_align(:,:,:,:,2));
% %--- Clear the last entry
% fprintf(repmat('\b', 1, lengthLastMsg));
% 
% %--- Display the progress of the sift-flow computation
% msg = '  Matrix decomposition (B-Channel)...';
% fprintf('%s',msg);
% lengthLastMsg = length(msg);
% pause(0.005);
% 
% [B_b,C_b] = matrix_decomposition(LR_LF_align(:,:,:,:,3));
% 
% %--- Clear the last entry
% fprintf(repmat('\b', 1, lengthLastMsg));
% 
% % Encode the principal basis to get the principal bases as an image
% [Ipb, param] = principal_basis_encoding(B_r(:,1), B_g(:,1), B_b(:,1),size(LR_LF_align));
% 
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % SUPERRESOLUTION OF THE PRINCIPAL BASIS
% %--------------------------------------------------------------------------
% 
% addpath('MATLAB\matconvnet\');
% addpath('MATLAB\matconvnet\matlab\');
% addpath('MATLAB\superresolution\srcnn\');
% 
% %run matconvnet/matlab/vl_setupnn;
% run vl_setupnn
% 
% model_filename = sprintf('MATLAB\\superresolution\\srcnn\\mf_%d_model.mat',mf);
% load(model_filename);
% 
% % Test that the convolution function works properly
% try
%     vl_nnconv(single(1),single(1),[]) ;
% catch
%     warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
%     vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose, ...
%         'enableImreadJpeg', false) ;
% end
% 
% % Super-resolve the current sub-aperture image
% Ipb_sr = srcnn(uint8(Ipb*255),net);
% 
% % Convert the image in range [0,1]
% Ipb_sr = double(Ipb_sr)/255;
% 
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % LIGHT FIELD RECONSTRUCTION
% %--------------------------------------------------------------------------
% 
% % Decode the principal basis
% [B_r(:,1), B_g(:,1),B_b(:,1)] = principal_basis_decoding(Ipb_sr,param);
% 
% % Initialize the SR_LF_align matrix
% SR_LF_align = zeros(size(LR_LF_align));
% 
% % Reconstruct the super-ersolved aligned light field
% SR_LF_align(:,:,:,:,1) = reshape(reshape(B_r*C_r,[size(LR_LF_align,1), ...
%     size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
%     [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
%     size(LR_LF_align,4)]);
% SR_LF_align(:,:,:,:,2) = reshape(reshape(B_g*C_g,[size(LR_LF_align,1), ...
%     size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
%     [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
%     size(LR_LF_align,4)]);
% SR_LF_align(:,:,:,:,3) = reshape(reshape(B_b*C_b,[size(LR_LF_align,1), ...
%     size(LR_LF_align,2),size(LR_LF_align,3)*size(LR_LF_align,4)]), ...
%     [size(LR_LF_align,1),size(LR_LF_align,2),size(LR_LF_align,3), ...
%     size(LR_LF_align,4)]);
% 
% % % Permute the channels of the aligned light field
% SR_LF_align = uint8(permute(SR_LF_align,[1,2,5,3,4]));

%--------------------------------------------------------------------------

% %--------------------------------------------------------------------------
% % LIGHT FIELD RECONSTRUCTION
% %--------------------------------------------------------------------------
% % Restore the original disparities using forward warping
% SR_LF = forward_warping(SR_LF_align,u,v);

clc; close all; clear all;
load temp

SR_LF = light_field_inpainting(SR_LF);

function X = light_field_inpainting(X_)

% Derive the mask that needs to be inpainted
mask = reshape(X_(:,:,1,:,:),[size(X_,1),size(X_,2),size(X_,4),size(X_,5)]) == -1; 

% Convert X_ to uint8
X = uint8(X_);

% Determine the number of missing pixels in each view
Nholes = reshape(sum(sum(mask,1),2),[size(X,4),size(X,5)]);

% Clear variables
clearvars X_;

% Derive the index of the smallest number of missing pixels
[~,idx] = sort(reshape(Nholes,[size(Nholes,1)*size(Nholes,2),1]),'ascend');

% Convert the index to 2D indices
[c_idx,r_idx] = ind2sub(size(Nholes),idx);

for i = 1:size(r_idx,1)
    if Nholes(r_idx(i),c_idx(i)) > 0 % This sub-aperture image must be restored
        % Compute the distance between the current sub-aperture image and
        % those which do not have holes
        dist = sqrt((c_idx(1:i-1) - c_idx(i)).^2 + (r_idx(1:i-1) - r_idx(i)).^2);
        
        % Find the minimum distance
        min_val = min(dist);
        
        % Find the indices of the minimum distances to be used as
        % references
        Jidx = find(dist == min_val);
        
        % Initialize the reference
        Xref = zeros(size(X,1),size(X,2),size(X,3),size(Jidx,1));
        
        % Extract the reference views
        for j = 1:size(Jidx,1)
            % Derive the reference images
            Xref(:,:,:,j) = X(:,:,:,r_idx(Jidx(j)),c_idx(Jidx(j)));
        end
        
        % Extract the sub-aperture image being processed
        I = X(:,:,:,r_idx(i),c_idx(i));
        
        % Extract the mask
        mk = mask(:,:,r_idx(i),c_idx(i));
        
        % Restore the current sub-aperture image
        I = inpainting(I, Xref,mk);
        
        % Put the restored sub-aperture image in the light field
        X(:,:,:,r_idx(i),c_idx(i)) = I;
    end
end

function Y = inpainting(I, Xref, mask)



