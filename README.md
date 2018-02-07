# Light Field Editing

The Light Field Editing repository was created to include scripts of a number of light field editing algorithms. This repository will include code covering different light field editing tools including light field spatial super-resolution, angular super-resolution, light field inpainting and light field denoising. This repository will include recent algorithms found in literature including methods being developed by our team.

Light Field Super-Resolution
----------------------------

The LF_super_resolution_analysis script is a command line script that allows to analyse the performance of a number of light field
spatial super-resolution algorithms. This script considers that the light field is down-sampled by a scale factor defined by mf and up-scaled to the target resolution using bi-cubic interpolation. The syntax to call this from the command line is:

LF_super_resolution_analysis('bicubic',3,true

Input: 

sr_method: specifies the super-resolution algorithm to be simulated. The following is a list of sr_methods that are supported here:

	- 'bicubic':	classical bicubic interpolation of each sub-				aperture image independently

	- 'lf_srcnn':	this method applies SRCNN to restore each 				sub-ape[rture image separately from the 					others [1],[2]
 
	- 'pca_rr': 	this method applies the pca_rr published
                        in [4] - super-resolves patch volumes
	
	- 'bm_pca_rr':	this method applies the pm_pca_rr published 				in [4] - super-resolves aligned patch 					volumes.

	- 'pb-srcnn': 	this method is the proposed method using
                        SRCNN to restore the principal basis

	- 'pb-vdsr': 	this method is the proposed method using VDSR 				to restore the principal basis
	- 'pb-lab402';    this mehtod is the proposed method using 					lab402 method to restore the principal basis

mf: numeric value that stands for the magnification factor that the method has to super-resolve

out_flag: This is a boolean value which specifies if the result attined in the simulation will be stored or not. By default it is set to false

[1] Y. Yoon, H. G. Jeon, D. Yoo, J. Y. Lee and I. S. Kweon, "Light-Field Image Super-Resolution Using Convolutional Neural Network," in IEEE Signal Processing Letters, vol. 24, no. 6, pp. 848-852, June 2017.
           
[2] C. Dong, C. C. Loy, K. He and X. Tang, "Image Super-Resolution Using Deep Convolutional Networks," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 38, no. 2, pp. 295-307, Feb. 1 2016.

[3] J. Kim, J. K. Lee and K. M. Lee, "Accurate Image Super-Resolution Using Very Deep Convolutional Networks," 2016 IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Las Vegas, NV, 2016, pp. 1646-1654.

[4] R.A. Farrugia, C. Galea, C. Guillemot, "Super Resolution of Light Field Images using Linear Subspace Projection of Patch-Volumes," in IEEE Journal on Selected Topics in Signal Processing, vol. 11, no. 7, pp. 1058-1071, Oct. 2017

Installation
------------

* You have to install matconvnet and include the matconvnet library into the folder MATLAB/
* Store all the light fields in a directory ../LF-DATASET/. We consider here four different datasets
	- EPFL (https://jpeg.org/plenodb/lf/epfl/)
	- INRIA (https://www.irisa.fr/temics/demos/IllumDatasetLF/index.html)
	- HCI (http://hci-lightfield.iwr.uni-heidelberg.de/)
	- STANFORD (http://lightfield.stanford.edu/lfs.html)
