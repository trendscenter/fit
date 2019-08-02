function [status, message] = ica_fuse_errorCheck(varIn, varCriteria)
% do error check for the variable with the mentioned criteria and
% return status and message.
%
% status 1 means error check went well
% status 0 means error occurred and see message variable for that

status = 1; message = '';

% define error criteria below
switch lower(varCriteria)

    case 'character'
        % check for valid character
        try
            [varOut] = strread(varIn, '%c');
        catch
            status = 0;
            message = 'Not a valid character';
        end

    case 'integer'
        % check if the variable is integer or not

        try
            varOut = strread(varIn, '%d');
        catch
            status = 0;
            message = 'Not a valid integer';
        end

    case 'double'
        % test whether the number entered is a valid double
        try
            varOut = strread(varIn, '%f');
        catch
            status = 0;
            message = 'Not a valid number';
        end

    case 'cell'
        if ~iscell(varIn)
            status = 0;
            error('Not a valid cell array');
        end

    case 'struct'
        if ~isstruct(varIn)
            status = 0;
            error('Not a valid structure');
        end

    case 'output_prefix'
        % check for valid output prefix

        if ~isempty(varIn)
            %[varOut] = ica_fuse_check_char(varIn);
            %             if isempty(varOut)
            %                 status = 0;
            %                 message = ['Not a valid output prefix. Don''t use characters like \, /, ?, :, ", <, > in the output prefix'];
            %             end
            try
                setappdata(0, varIn, 1);
                rmappdata(0, varIn);
                status = 1;
            catch
                status = 0;
                message = ['Specified output prefix (', varIn, ') will not generate a valid name to store the application data'];
            end
        end

    case 'file'
        % check if the file exists or not
        if ~exist(varIn, 'file')
            status = 0;
            message = ['file: ', varIn , ' doesn''t exist'];
        end

    case 'directory'
        % check if the directory exists or not
        if ~exist(varIn, 'dir')
            status = 0;
            message = ['directory: ', varIn , ' doesn''t exist'];
        end

    otherwise
        status = 0;
        message = 'Unknown error check criteria for the variable specified.';
end