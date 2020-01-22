function inputVector = ica_fuse_unitRamp(inputVector)

% calculates the unit ramp function
maxVal = max(inputVector);

inputVector = inputVector/maxVal;