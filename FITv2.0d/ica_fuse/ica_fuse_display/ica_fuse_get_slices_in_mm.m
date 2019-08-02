function [slices_in_mm, displayParameters] = ica_fuse_get_slices_in_mm(displayParameters, anatomicalView)

if ~isfield(displayParameters, 'structVol')
    ica_fuse_defaults;
    global ANATOMICAL_FILE;
    structVol = ica_fuse_getVol(ANATOMICAL_FILE, 1);
    displayParameters.structVol = structVol;
end

[sliceParameters] = ica_fuse_get_slice_def(displayParameters.structVol, anatomicalView);
slices_in_mm = sliceParameters.slices; clear sliceParameters;
slices_in_mm = ica_fuse_constructString(slices_in_mm);