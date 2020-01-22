function compData = ica_fuse_applyDispParameters(plotType, compData, convertToZ, returnValue, threshValue)
% apply display parameters to the component images

threshValue = abs(threshValue);

if strcmpi(plotType, 'image')
    
    [size_data] = size(compData);
    
    % reshape component image
    compData = reshape(compData, size_data(1), prod(size_data(2:end)));
    
    % convertToZScores
    if convertToZ
        compData = ica_fuse_convertImageToZScores(compData);
    end
    
    %get desired image values
    if returnValue == 1
        compData = compData;
    elseif returnValue == 2
        temp = zeros(size(compData));
        indices = find(compData>0);
        temp(indices) = compData(indices);
        compData = temp;
        clear temp;
    elseif returnValue == 3
        compData = abs(compData);
        
        % negative image values
    elseif returnValue == 4
        temp = zeros(size(compData));
        indices = find(compData < 0);
        temp(indices) = compData(indices);
        compData = temp;
        clear temp;
    end
    
    %threshold image
    for i=1:size(compData,1)
        compData(i, :) = applyThreshValue(compData(i, :), returnValue, threshValue);
    end
    
    % reshape the image to 3D
    compData = reshape(compData, size_data);
    
end



function data = applyThreshValue(data, returnValue, threshold)

if (length(threshold) > 1)
    
    if ((returnValue == 2) || (returnValue == 3))
        % Positive or absolute
        data(data < min(threshold)) = 0;
        data(data > max(threshold)) = max(threshold);
    elseif (returnValue == 4)
        % Negative
        data(abs(data) < min(threshold)) = 0;
        data(abs(data) > max(threshold)) = -max(threshold);
    else
        % Positive and Negative
        neg_threshold = -threshold;
        
        pos_inds = (data > 0);
        neg_inds = (data < 0);
        
        % Handle positive range
        tmp1 = data(pos_inds);
        tmp1(tmp1 < min(threshold)) = 0;
        tmp1(tmp1 > max(threshold)) = max(threshold);
        
        % handle negative range
        tmp2 = data(neg_inds);
        tmp2(tmp2 > max(neg_threshold)) = 0;
        tmp2(tmp2 < min(neg_threshold)) = min(neg_threshold);
        
        data(pos_inds) = tmp1;
        data(neg_inds) = tmp2;
    end
    
    %data(abs(data) < min(threshold)) = 0;
    %data(abs(data) > max(threshold)) = 0;
else
    data(abs(data) < threshold) = 0;
end