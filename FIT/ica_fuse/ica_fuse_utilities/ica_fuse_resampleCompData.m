function [compData] = ica_fuse_resampleCompData(compData, interpolatedLength, voxels, dataLength)
% Resample component data using resample function.
%
% Input:
% 1. compData - data to be sampled
% 2. interpolatedLength - Interpolated data length
% 3. voxels - factor to interpolate
% 4. dataLength - factor to decimate
%
% Output:
% compData - resampled data


interpFactor = ceil(voxels/dataLength);

if interpFactor > 1

    if size(compData, 1) ~= interpolatedLength
        temp = zeros(interpolatedLength, 2, size(compData, 3));

        % Loop over number of subjects
        for ii = 1:size(temp, 3)
            % Resample data
            temp(:, 2, ii) = ica_fuse_resample(compData(:, 2, ii), voxels, dataLength);
        end
        % end loop over number of subjects

        % Resample x axis
        xAxis = ica_fuse_resample(compData(:, 1, 1), voxels, dataLength);
        for ii = 1:size(temp, 3)
            % Resample data
            temp(:, 1, ii) = xAxis;
        end
        clear compData;
        compData = temp;
        clear temp;
        clear xAxis;
    end
    
end