function [varargout] = ica_fuse_downsample(data, timeAxis, downSampleFactor)
% Downsample data and time axis using downsample function (Signal
% Processing Toolbox)
%
% Input:
% 1. data - data to be downsampled
% 2. timeAxis - timeAxis to be downsampled
% 3. downSampleFactor - factor to downsample
%
% Output:
% 1. data - Downsampled data
% 2. timeAxis - Time axis

if ~exist('data', 'var')
    error('Data variable must be passed');
end

if ~exist('timeAxis', 'var')
    timeAxis = [];
end


temp1 = data;
temp2 = timeAxis;

checkFunc = which('downsample.m');

if ~isempty(checkFunc)
    try
        % Data
        data = downsample(data, downSampleFactor);
        if ~isempty(timeAxis)
            % Time axis
            timeAxis = downsample(timeAxis, downSampleFactor);
        end
    catch
        data = temp1;
        timeAxis = temp2;
    end
end

% Return output
if isempty(timeAxis)
    varargout{1} = data;
else
    varargout{1} = data;
    varargout{2} = timeAxis;
end
