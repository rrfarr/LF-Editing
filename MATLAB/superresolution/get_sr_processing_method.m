function proc_method = get_sr_processing_method(sr_method)
% This function maps the processing method from the super-resolution method
% considered. We are considering here two possibilities:
%
% 1 - independentSR - super-resolution is applied on each sub-aperture
% image independent from the other views. Any single image super-resolution
% or interpolation can be used here
%
% 2 - coherentSR - the light field is considered as a volume and
% super-resolution is performed by exploiting the light field structure.
%

if strcmp(sr_method, 'bicubic')
    proc_method = 'independentSR';
elseif strcmp(sr_method,'srcnn') || strcmp(sr_method,'lf_srcnn')
    proc_method = 'independentSR';
elseif strcmp(sr_method,'vdsr')
    proc_method = 'independentSR';
end

