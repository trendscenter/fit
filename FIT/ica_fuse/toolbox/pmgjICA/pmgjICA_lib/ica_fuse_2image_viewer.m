function [cData, ret_icasig] = ica_fuse_image_viewer(file_names, savepath, varargin)


%%% Example
%%im_data1=ica_fuse_image_viewer(files, savepath, 'labels', labels, 'display_type', 'montage', 'structfile', ...
%%        fullfile(fileparts(which('gift.m')), 'ica_fuse_templates', 'ch2bet.nii'), ...
%%        'threshold', thr, 'slices_in_mm', (-40:4:72), 'convert_to_zscores', ...
%%        convert_to_zscores, 'image_values', 'positive and negative','iscomposite','no', 'savefig', 'yes');

%%% old function arguments
%%%%%function dispParameters= ica_fuse_image_viewer(file_names, varargin);   
%% Image viewer
%
% Inputs:
% 1. file_names - Character array of file names
% Arguments passed in pairs
%   a. display_type - montage, render or ortho
%   b. Image_values - Positive and Negative, Positive, Absolute, Negative
%   c. anatomical view - Axial, sagittal or coronal (for montage)
%   d. convert_to_zscores - Convert to z-scores
%   e. threshold - Threshold, including both a lower and an optional upper limit (e.g., 2 if only a lower limit or [2, 8], having both lower and upper limit)
%
%  example command, where file may be a string to a nii with several spatial components: ica_fuse_image_viewer(file, 'display_type', 'montage', 'structfile',fullfile(fileparts(which('gift.m')), 'ica_fuse_templates', 'ch2bet.nii'), 'threshold', [2, 8], 'slices_in_mm', (-40:4:72), 'convert_to_zscores','yes', 'image_values', 'positive and negative','yes','iscomposite','yes');

%% use save id, even if it not available you have to use as 0 or 1, same as you may use save_path as empty it is not used

image_values = 'positive and negative';
threshold = 0;
convert_to_zscores = 'no';
structFile = [];
display_type = 'montage';
anatomical_view = 'axial';
slices_in_mm = [];
isComposite = 'no';
colorbar_title = '';
FONT_COLOR = [1 1 1]; %ce042525


ica_fuse_defaults;
global UI_FS;
global UI_FONTNAME;
global TEXT_DISPLAY_SLICES_IN_MM;
global FONT_COLOR;

labels = '';
useColorbar = 1;

%% Loop over variable no. of arguments
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'structfile'))
        structFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'image_values'))
        image_values = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'convert_to_zscores'))
        convert_to_zscores = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'labels'))
        labels = varargin{n + 1};
        %     elseif (strcmpi(varargin{n}, 'fig_title'))
        %         fig_title = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'display_type'))
        display_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'slices_in_mm'))
        slices_in_mm = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'anatomical_view'))
        anatomical_view = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'axesh'))
        axesH = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar'))
        useColorbar = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'iscomposite'))
        isComposite = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'usegui'))
        useGUI = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'colorbar_title'))
        colorbar_title = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'savefig'))
        savefig = varargin{n+1};
    elseif(strcmpi(varargin{n}, 'save_nifti'))
        save_nifti = varargin{n+1};
    end
end

if (~exist('useGUI', 'var'))
    useGUI = 0;
    if (~exist('file_names', 'var'))
        useGUI = 1;
    end
end

if (useGUI)
    
    if (~exist('file_names', 'var'))
        file_names = ica_fuse_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.img;*.nii', 'filetype', 'image', 'title', ...
            'Select file/files ...');
    end
    
    drawnow;
    
    if (isempty(file_names))
        error('File is not selected');
    end
    
    getDispDefaults;
    
    appStr = 'disp_para_data';
    if (isappdata(0, appStr))
        dispParameters = getappdata(0, appStr);
        rmappdata(0, appStr);
    else
        error('Input parameters window was quit');
    end
    
    
    structFile = dispParameters.structFile;
    display_type = dispParameters.display_type;
    convert_to_zscores = dispParameters.convert_to_z;
    threshold = dispParameters.threshold;
    image_values = dispParameters.image_values;
    anatomical_view = dispParameters.slice_plane;
    slices_in_mm = dispParameters.slices_in_mm;
    isComposite = dispParameters.isComposite;
    
end

if (isnumeric(isComposite))
    compositeOpts = {'No', 'Yes'};
    isComposite = lower(compositeOpts{isComposite});
end


file_names = char(file_names);

% Call function recursively
if (strcmpi(isComposite, 'yes'))
    ica_fuse_display_composite(file_names, 'anatomical_file', structFile, 'image_values', image_values, 'threshold', threshold, 'convert_to_zscores', convert_to_zscores, 'labels', labels, ...
        'display_type', display_type, 'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_view, 'colorbar', useColorbar);
    return;
end
%     for nF = 1:size(file_names, 1)
%         ica_fuse_image_viewer(squeeze(file_names(nF, :)), 'structfile', structFile, 'image_values', image_values, 'threshold', threshold, 'convert_to_zscores', convert_to_zscores, 'labels', labels, ...
%             'display_type', display_type, 'slices_in_mm', slices_in_mm, 'anatomical_view', anatomical_view, 'colorbar', useColorbar);
%     end

file_names = ica_fuse_rename_4d_file(file_names);
%file_names = deblank(file_names(1, :));

if (isempty(structFile))
    structFile = deblank(file_names(1, :));
end

allFileNames = file_names;

for nF = 1:size(allFileNames, 1)
    
    file_names = deblank(allFileNames(nF, :));
    
    fontSizeText = 8;
    convertToZ = strcmpi(convert_to_zscores, 'yes');
    returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');    
    
    if (iscell(labels) && file_names(end:end)~='1')
        currentLabel = labels{nF};
        fig_title = currentLabel;
    elseif (iscell(labels) && file_names(end:end)=='1')
        disp('Inside......')
        currentLabel = labels;
        fig_title = currentLabel{nF};
    else
        [dd, FileName, extn] = fileparts(deblank(file_names));
        [dd, label_number] = ica_fuse_parseExtn(extn);
        currentLabel = [FileName, ',', num2str(label_number)];
        fig_title = currentLabel;
    end
    
    
    %     if (isempty(labels))
    %         [dd, FileName, extn] = fileparts(deblank(file_names(1, :)));
    %         [dd, label_number] = ica_fuse_parseExtn(extn);
    %         labels = [FileName, ',', num2str(label_number)];
    %     end   
    
    
    %load ica_fuse_colors coldhot_sensitive;
    load ica_fuse_colors coldhot cold hot;
    if (returnValue == 1)
        % positive and negative
        cmap = coldhot(1:4:end, :);
    elseif (returnValue == 4)
        % negative
        cmap = cold(1:4:end, :);
    else
        % hot
        cmap = hot(1:4:end, :);
    end
    %     if (returnValue == 1)
    %         cmap = coldhot_sensitive(1:4:end, :);
    %     elseif (returnValue == 4)
    %         cmap = coldhot_sensitive(1:128, :);
    %         cmap = cmap(1:2:end, :);
    %     else
    %         cmap = coldhot_sensitive(129:end, :);
    %         cmap = cmap(1:2:end, :);
    %     end
    
    
    if (strcmpi(display_type, 'montage'))
        %% Montage
        
        if (~exist('axesH', 'var'))
            
            gH = ica_fuse_getGraphics(fig_title, 'graphics', 'Image Viewer', 'on');
            axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        end
        
        structVol = ica_fuse_spm_vol(structFile);
        
        if (isempty(slices_in_mm))
            sP = ica_fuse_get_slice_def(structVol, anatomical_view);
            slices_in_mm = sP.slices;
        end
        
        % resize the image and return header Info
        [images, coords, HInfo, slices_in_mm, text_left_right] = ica_fuse_resizeImage(structVol, file_names, anatomical_view, ...
            slices_in_mm, 1);
        %%disp('images...')
        %%size(images)  %%% 2   182   218    29
        % get the structural image
        structuralImage = squeeze(images(1, :, :, :));
        
        % get the component images
        icasig = images(2:end, :, :, :);        
        size(icasig)

        clear images;
        
        % [structuralImage, icasig, coords, HInfo] = ica_fuse_returnResizedImage(structVol, file_names, anatomical_view, slices_in_mm, 'real', [], 1);
        structDIM = HInfo.DIM
        icasig = reshape(icasig, 1, prod(structDIM));
        % apply display parameters
        icasig = ica_fuse_2applyDispParameters(icasig, convertToZ, returnValue, threshold, structDIM, HInfo);     
        
        ret_icasig=icasig;  %%% 182   218    29
        %%ret_icasig = reshape(ret_icasig, structDIM);  %%% %%% 182   218    29   %%% before montage
        %%save('icasig.mat', 'ret_icasig')

        size(ret_icasig)


        %get colormap
        cmap = [cmap; gray(size(cmap, 1))];
        %cmap = ica_fuse_getColormap(1, returnValue, 1);
        
        % return overlayed images depending upon the data type of icasig and
        % structuralImage
%         [icasig,  maxICAIM, minICAIM, minInterval, maxInterval] = ica_fuse_check_overlayComp(icasig, structuralImage, ...
%             returnValue, 1);

        icaDIM = [size(icasig, 1), size(icasig, 2), size(icasig, 3)];
        

        %%%% Note if we provide range  int threshold, like [1, 6], range
        %%%% min/max 6. then maxICAIM and minICAIM will 6/-6 if the voxel
        %%%% values lower than the range.

        [icasig,  maxICAIM, minICAIM, minInterval, maxInterval] = ica_fuse_2overlayImages(icasig,structuralImage, icaDIM, structDIM, returnValue, threshold);       
        icasig = reshape(icasig, structDIM);
        class(icasig)        
        size(icasig)
        structVol
        if strcmpi(save_nifti, 'yes')
            disp('Nothing here ...')
            %%save('icasig.mat', 'icasig')
            
            %%info=niftiinfo('smtest22.nii');
            %%info.ImageSize=structDIM;
            %%info.Filename=savepath{nF};
            %%niftiwrite(icasig, savepath{nF}, info);

            % % %             A=ica_fuse_spm_vol('smtest22.nii');
            % % %             A.fname=savepath{nF};
            % % %             A.dim=structDIM;
            % % %             A=rmfield(A,'pinfo');
            % % %             A=rmfield(A,'descrip');
            % % %             A=rmfield(A,'dt');
            % % %             A=rmfield(A,'n');
            % % %             ica_fuse_write_nifti_data(savepath{nF}, A, ret_icasig);

            %%A=ica_fuse_spm_vol(structFile);
            %%%%%data=ica_fuse_loadData('C:\Users\oagcaoglu\Downloads\icasig.mat');
            % %             A=ica_fuse_spm_vol(fullfile(fileparts(which('gift.m')), 'ica_fuse_templates', 'ch2bet.nii'));
            % %             %%disp('Inside save_niftii........')
            % %             %%%size(ret_icasig)
            % %             save('mytmp.mat', 'ret_icasig');
            % %             ret_icasig=ica_fuse_loadData('mytmp.mat');
            % %             %%%size(ret_icasig)
            % % 
            % %             data_comb=reshape(ret_icasig,182,218,1,29);
            % %             size(data_comb)
            % %             A.fname=savepath{nF};
            % %             A=rmfield(A,'pinfo');
            % %             A=rmfield(A,'descrip');
            % %             A=rmfield(A,'dt');
            % %             A=rmfield(A,'n');
            % %             sz=size(data_comb);
            % %             A.dim = sz(1:3);
            % %             for i=1:size(data_comb,4)
            % %                 VV(i,1)=A;
            % %             end
            % %             ica_fuse_write_nifti_data(A.fname, VV, data_comb);
            % % 
            % %             clear VV
            % %             
            % %             %%%reread=ica_fuse_loadData('my_icasig.nii');
            % % 
            % %             %%% checking after load
            % %             %%a=abs(reread-data_comb);
            % %             %%%sum(a(:))
            % %             %%%max(a(:))
            % %             save_png=strrep(savepath{nF},'.nii','.png');
            % %             save_create_2d(ret_icasig, save_png);

        end
        
        if strcmpi(anatomical_view, 'sagittal')
            icasig = permute(icasig, [2, 3, 1]);
            HInfo.VOX = [HInfo.VOX(2) HInfo.VOX(3) HInfo.VOX(1)];
        end
        
        if strcmpi(anatomical_view, 'coronal')
            icasig = permute(icasig, [1, 3, 2]);
            HInfo.VOX = [HInfo.VOX(1) HInfo.VOX(3) HInfo.VOX(2)];
        end
        
        DIM(1) = size(icasig, 1);
        DIM(2) = size(icasig, 2);
        DIM(3) = size(icasig, 3);
        
        
        
        [im, numImagesX, numImagesY, slices_in_mm_new] = ica_fuse_2returnMontage(icasig, [], DIM, HInfo.VOX, slices_in_mm);
        
        
        colormap(cmap);
        
        ImageAxis = image(im, 'parent', axesH, 'CDataMapping', 'scaled');
        set(axesH, 'clim',  [minInterval, 2*maxInterval]); % set the axis positions to the specified
        axis(axesH, 'off');
        axis(axesH, 'image');
        
        
        %if (exist('labels', 'var'))
        %% original code
        %%title(currentLabel, 'parent', axesH, 'color', FONT_COLOR, 'Horizontalalignment', 'center'); 
        %end
        if (iscell(labels) && file_names(end:end)=='1')
            FONT_COLOR = [1 1 1]; %ce042525
            title(currentLabel, 'color', FONT_COLOR, 'FontSize', 15);
        else
            FONT_COLOR = [1 1 1]; %ce042525
            title(currentLabel, 'parent', axesH, 'color', FONT_COLOR, 'Horizontalalignment', 'center'); 
        end
        
        if strcmpi(TEXT_DISPLAY_SLICES_IN_MM, 'on')
            dim(1) = size(im, 1);
            dim(2) = size(im, 2);
            % name the text and place it in the correct order (slices in mm).
            textCount = 0;
            yPos = 1 + dim(1) / numImagesY;
            for nTextRows = 1:numImagesY
                xPos = 1;
                for nTextCols = 1:numImagesX
                    textCount = textCount + 1;
                    if textCount <= DIM(3)
                        txHandle(textCount) = text(xPos, yPos, num2str(round(slices_in_mm_new(textCount))), 'color', FONT_COLOR,  ...
                            'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                            'FontName', UI_FONTNAME, 'parent', axesH);
                    end
                    xPos = xPos + (dim(2) / numImagesX);
                end
                % end for cols
                yPos = yPos + (dim(1) / numImagesY); % update the y position
            end
            % end for rows
            
        end
        
        if (useColorbar)
            
            imagePos = get(axesH, 'position');
            pos = [imagePos(1) + imagePos(3) + .04, imagePos(2)+.1, 0.03, 0.13];
            
            if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
                ColorbarHandle = colorbar('peer', axesH);
            else
                ColorbarHandle = colorbar;
            end
            set(ColorbarHandle, 'position', [pos(1), pos(2) pos(3) pos(4)]);
            
            
            maxLabel = round(maxICAIM(1)*10)/10;
            minLabel = round(minICAIM(1)*10)/10;

         

            
            maxLabel = char(sprintf('%0.1f', maxLabel));
            minLabel = char(sprintf('%0.1f', minLabel));
            
            minLabel= padStrs(minLabel, colorbar_title);
            
            YTick = [minInterval, maxInterval];
            set(ColorbarHandle,'YLim',YTick);
            
            set(ColorbarHandle, 'YTick', []);
            
            if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
                
                title(maxLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
                xlabel(minLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
                
            else
                
                labelsH = get(ColorbarHandle, 'label');
                set(ColorbarHandle, 'color', FONT_COLOR);
                set(labelsH, 'units', 'normalized', 'fontsize', 15);
                title(maxLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
                set(labelsH, 'string', minLabel, 'rotation',0, 'fontsize', 15);
                set(labelsH, 'position', [0.45, -.12, 0]);
                
            end
            
            set(ColorbarHandle, 'YTickLabel', []);
            
        end
        
    elseif (strcmpi(display_type, 'ortho'))
        %% Ortho Plot
        
        gH = ica_fuse_orth_views(file_names, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_zscores, 'threshold', threshold, 'set_to_max_voxel', 1, ...
            'labels', '', 'fig_title', fig_title, 'cmap', cmap);
        
    elseif (strcmpi(display_type, 'render'))
        %% Render
        
        [V, HInfo] = ica_fuse_returnHInfo(file_names);
        [R, C, P] = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
        RCP = [R(:)';C(:)';P(:)'];
        data = ica_fuse_spm_read_vols(V);
        data(isfinite(data) == 0) = 0;
        data = ica_fuse_2applyDispParameters(data(:)', convertToZ, returnValue, threshold, V.dim(1:3), HInfo);
        mask = find(abs(data) > eps);
        
        XYZ = RCP(:, mask);
        Z = data(mask); %ones(length(mask), 1);
        dat   = struct(	'XYZ',	XYZ,...
            't',	Z, ...
            'mat',	V.mat, ...
            'dim',	V.dim(1:3)');
        global prevrend;
        rendfile = fullfile(fileparts(which('groupica.m')), 'ica_fuse_templates', 'render_single_subj.mat');
        prevrend.rendfile = rendfile;
        prevrend.col = [1 0 0; 0 1 0; 0 0 1];
        prevrend.brt = NaN;
        [rend, ddd, cmap] = ica_fuse_spm_render(dat, NaN, rendfile);
        rend = concatDat(rend);
        
        
        if (~exist('axesH', 'var'))
            gH = ica_fuse_getGraphics(fig_title, 'graphics', 'Image Viewer', 'on');
            axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        end
        
        ImageAxis = image(rend, 'parent', axesH);
        axis(axesH, 'image');
        axis(axesH, 'off');
        
        cmap(1:64, :) = [];
        
        %% Put colorbar
        if (~isempty(cmap))
            
            colormap(cmap);
            ch = colorbar('peer', axesH);
            pos = get(ch, 'position');
            pos(4) = 0.13;
            pos(2) = 0.5 - 0.5*pos(4);
            pos(1) = pos(1) + 0.05;
            pos(3) = 0.03;
            set(ch, 'position', pos);
            set(ch, 'YTick', []);
            
            minVal = min(Z(:));
            maxVal = max(Z(:));
            
            if (minVal < 0)
                if (maxVal > 0)
                    minVal = -maxVal;
                end
            end
            
            minLabel = char(sprintf('%0.1f', minVal));
            maxLabel = char(sprintf('%0.1f', maxVal));
            
            if ((ica_fuse_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
                title(maxLabel, 'parent', ch, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
                xlabel(minLabel, 'parent', ch, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
            else
                labelsH = get(ch, 'label');
                set(ch, 'color', FONT_COLOR, 'fontsize', 15);
                set(labelsH, 'units', 'normalized', 'fontsize', 15);
                title(maxLabel, 'parent', ch, 'color', FONT_COLOR, 'Horizontalalignment', 'center', 'fontsize', 15);
                set(labelsH, 'string', minLabel, 'rotation',0);
                set(labelsH, 'position', [0.45, -.12, 0]);
            end
        end
        
    end

    if (strcmpi(savefig, 'yes'))

       
        print(gcf,savepath{nF}, '-dpng'); %% OK

        %saveas(gcf, savepath{nF}); %% save for png, %%% OK

        cData=print('-RGBImage');
        %%%all_figs(:,:,:, nF)=CData;

        %%%% montage in pdf
        %%%figg=gcf;
        %%%%figg
        %%%%allimages(:,:,:,nF)=print('-RGBImage');
        %%%%fig=gcf;
        %%%fig = get(groot,'CurrentFigure');

        %exportgraphics(axesH,sfilename,'BackgroundColor','current','Append',true);
        % %%%% OK
        %close(gcf);
    end
    
    clear axesH;
    
end

function MAT = concatDat(dat)
%% Concatenate rendered data
%

nCols = 2;
nLoops = ceil(length(dat)/nCols);

row_size = 0;
col_size = [];
for nLoop = 1:nLoops
    dims = size(dat{2*nLoop-1});
    row_size = row_size + dims(1);
    col_size = max([col_size, dims(2)]);
end

col_size = 2*col_size;

MAT = zeros(row_size, col_size, 3);

er = 0;
loopNum = 0;

for i = 1:nLoops
    
    loopNum = loopNum + 1;
    ec = 0;
    if (loopNum <= length(dat))
        tmp = dat{loopNum};
        
        current_dims = size(tmp);
        
        sr = er + 1;
        er = er + current_dims(1);
        
        sc = ec + round((col_size/2 - current_dims(2))/2)+1;
        ec = sc + current_dims(2)-1;
        
        MAT(sr:er, sc:ec, :) = tmp;
    end
    
    loopNum = loopNum + 1;
    
    if (loopNum <= length(dat))
        tmp = dat{loopNum};
        
        current_dims = size(tmp);
        
        sc = round(col_size/2 + (col_size/2 - current_dims(2))/2)+1;
        ec = sc + current_dims(2)-1;
        
        MAT(sr:er, sc:ec, :) = tmp;
    end
    
end



function  getDispDefaults
%% Get display defaults
%

ica_fuse_defaults;
global UI_FS;
global UI_FONTNAME;
global CONVERT_Z;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;

%% Draw graphics
figureTag = 'setup_image_viewer_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = ica_fuse_getGraphics('Image viewer Params', 'displaygui', figureTag, 'on');
set(InputHandle, 'menubar', 'none');

yPos = 0.92;
promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight;
okWidth = 0.12; okHeight = promptHeight;

%% Select display type
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select display type', 'tag', 'prompt_display_type', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', char('Montage', 'Ortho', 'Render'), 'tag', ...
    'display_type', 'fontsize', UI_FS - 1, 'value', 1, 'callback', {@displayTypeCallback, InputHandle});


%% Is it a composite plot?
compositeOptions = {'No', 'Yes'};
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Do You Want To Use Composite Plot?', 'tag', 'prompt_composite', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
chkInd = strmatch('no', lower(compositeOptions), 'exact');
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', {'No', 'Yes'}, 'tag', 'isComposite', 'fontsize', UI_FS - 1, ...
    'value', chkInd);

%% Z-scores
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Convert To Z-scores?', 'tag', 'prompt_z', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
zOptions = char('Yes', 'No');
chkInd = strmatch(lower(CONVERT_Z), lower(zOptions), 'exact');
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', zOptions, 'tag', 'convert_to_z', 'fontsize', UI_FS - 1, ...
    'value', chkInd);

%% Enter threshold
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter threshold', 'tag', 'prompt_threshold', 'fontsize', ...
    UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', THRESHOLD_VALUE, 'tag', 'threshold', 'fontsize', UI_FS - 1);


%% Image values
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Image Values', 'tag', 'prompt_display', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);


imOptions = char('Positive and Negative', 'Positive', 'Absolute value', 'Negative');
returnValue = strmatch(lower(IMAGE_VALUES), lower(imOptions), 'exact');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', imOptions, 'tag', 'image_values', 'fontsize', UI_FS - 1, ...
    'value', returnValue);


%% Slice plane
planeOptions = {'Axial', 'Sagittal', 'Coronal'};
planeVal = strmatch(lower(ANATOMICAL_PLANE), lower(planeOptions));

promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select anatomical plane', 'tag', 'prompt_slice_plane', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', planeOptions, 'tag', ...
    'slice_plane', 'fontsize', UI_FS - 1, 'value', planeVal, 'callback', {@updateSlicesInMM, InputHandle});


%% Slices in mm
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter slices in mm', 'tag', 'prompt_slice_range', ...
    'fontsize', UI_FS - 1);
ica_fuse_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', '', 'tag', ...
    'slices_in_mm', 'fontsize', UI_FS - 1, 'value', 1);


%% Load anatomical

startPath = fileparts(which('groupica.m'));
startPath = fullfile(startPath, 'ica_fuse_templates');
structFile = fullfile(startPath, 'ch2bet.nii');

anatWidth = 0.2;
anatHeight = 0.05;
anatPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, anatWidth, anatHeight];
anatPos(2) = anatPos(2) - 0.5*anatPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', anatPos, 'string', 'Select Anatomical', 'tag', 'anat_button', 'fontsize',...
    UI_FS - 1, 'callback', {@selectAnatomical, InputHandle}, 'userdata', structFile);


%% Add cancel, save and run buttons
okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
ica_fuse_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Ok', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@getDefs, InputHandle});

slicePlaneH = findobj(InputHandle, 'tag', 'slice_plane');

updateSlicesInMM(slicePlaneH, [], InputHandle);

waitfor(InputHandle);

function displayTypeCallback(hObject, event_data, figH)
%% Display type callback


displayTypeH = findobj(figH, 'tag', 'display_type');
displayTypeOptions = cellstr(get(displayTypeH, 'string'));
displayTypeOptionsVal = get(displayTypeH, 'value');
displayType = lower(deblank(displayTypeOptions{displayTypeOptionsVal}));

anatH = findobj(figH, 'tag', 'anat_button');

if (strcmpi(displayType, 'render'))
    set(anatH, 'enable', 'off');
else
    set(anatH, 'enable', 'on');
end


spH = findobj(figH, 'tag', 'slice_plane');

if (~strcmpi(displayType, 'montage'))
    set(spH, 'enable', 'inactive');
else
    set(spH, 'enable', 'on');
end

sliceInMMH = findobj(figH, 'tag', 'slices_in_mm');

if (~strcmpi(displayType, 'montage'))
    set(sliceInMMH, 'enable', 'inactive');
else
    set(sliceInMMH, 'enable', 'on');
end

function selectAnatomical(hObject, event_data, figH)
%% Anatomical callback
%

startPath = fileparts(which('groupica.m'));
startPath = fullfile(startPath, 'ica_fuse_templates');

oldDir = pwd;

if (~exist(startPath, 'dir'))
    startPath = pwd;
end

% get the structural file
structFile = ica_fuse_selectEntry('typeEntity', 'file', 'title', 'Select Structural File', 'filter', ...
    '*.img;*.nii', 'fileType', 'image', 'fileNumbers', 1, 'startpath', startPath);

drawnow;

cd(oldDir);

if (~isempty(structFile))
    set(hObject, 'userdata', structFile);
    sliceInMMH = findobj(figH, 'tag', 'slice_plane');
    updateSlicesInMM(sliceInMMH, [], figH);
end


function updateSlicesInMM(slicePlaneH, event_data, figH)
%% Update slices in mm
%

anatH = findobj(figH, 'tag', 'anat_button');
structFile = get(anatH, 'userdata');

% slice plane
%slicePlaneH = findobj(figH, 'tag', 'slice_plane');
sliceOptions = cellstr(get(slicePlaneH, 'string'));
sliceOptionsVal = get(slicePlaneH, 'value');
slicePlane = lower(deblank(sliceOptions{sliceOptionsVal}));

% % Display type
displayTypeH = findobj(figH, 'tag', 'display_type');
displayTypeOptions = cellstr(get(displayTypeH, 'string'));
displayTypeOptionsVal = get(displayTypeH, 'value');
displayType = lower(deblank(displayTypeOptions{displayTypeOptionsVal}));

drawnow;

% Compute slices in mm
imagVol = ica_fuse_2returnHInfo(structFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = ica_fuse_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = ica_fuse_constructString(slices_in_mm);

% slices in mm
sliceInMMH = findobj(figH, 'tag', 'slices_in_mm');
set(sliceInMMH, 'string', slices_in_mm);

if (~strcmpi(displayType, 'montage'))
    set(sliceInMMH, 'enable', 'inactive');
else
    set(sliceInMMH, 'enable', 'on');
end


function getDefs(hObject, event_data, handles)
%% get defs

% display type
displayH = findobj(handles, 'tag', 'display_type');
displayStr = cellstr(get(displayH, 'string'));
displayVal = get(displayH, 'value');
display_type = lower(displayStr{displayVal});
dispParameters.display_type = display_type;

% convert to z
convertToZH = findobj(handles, 'tag', 'convert_to_z');
zStr = cellstr(get(convertToZH, 'string'));
zVal = get(convertToZH, 'value');
convert_to_z = lower(zStr{zVal});
dispParameters.convert_to_z = convert_to_z;

% Threshold
dispParameters.threshold = str2num(get(findobj(handles, 'tag', 'threshold'), 'string'));

% Image values
imH = findobj(handles, 'tag', 'image_values');
imStr = cellstr(get(imH, 'string'));
imVal = get(imH, 'value');
image_values = lower(imStr{imVal});
dispParameters.image_values = image_values;

% Slice plane
spH = findobj(handles, 'tag', 'slice_plane');
spStr = cellstr(get(spH, 'string'));
spVal = get(spH, 'value');
slice_plane = lower(spStr{spVal});
dispParameters.slice_plane = slice_plane;

% Slices in mm
dispParameters.slices_in_mm = str2num(get(findobj(handles, 'tag', 'slices_in_mm'), 'string'));

% anatomical file
anatH = findobj(handles, 'tag', 'anat_button');
structFile = get(anatH, 'userdata');
dispParameters.structFile = structFile;

% is composite?
imH = findobj(handles, 'tag', 'isComposite');
imStr = cellstr(get(imH, 'string'));
imVal = get(imH, 'value');
isComposite = lower(imStr{imVal});
dispParameters.isComposite = isComposite;

setappdata(0, 'disp_para_data', dispParameters);

delete(handles);

drawnow;


function strs = padStrs(varargin)
% Pad Strings

strs = char(varargin{:});
strs = repmat(' ', size(strs));
number_chars = size(strs, 2);

for n = 1:length(varargin)
    current_char = varargin{n};
    midPoint = ceil(number_chars/2);
    startPoint = max([ceil(midPoint - 0.5*length(current_char)), 1]);
    strs(n, startPoint:startPoint + length(current_char)-1) = current_char;
end

function [minICAIM, maxICAIM] = set_ic_limits(minICAIM, maxICAIM, threshold, returnValue)
% Set IC limits
%

minICAIM = abs(minICAIM);
maxICAIM = abs(maxICAIM);
threshold = abs(threshold);

if ((returnValue == 2) || (returnValue == 3))
    % absolute or positive values
    minICAIM = min([minICAIM, threshold]);
    maxICAIM = max([maxICAIM, threshold]);
elseif (returnValue == 4)
    % negative values
    tmp_maxicaim = max([minICAIM, maxICAIM, threshold]);
    tmp_minicaim = min([minICAIM, maxICAIM, threshold]);
    minICAIM = -tmp_maxicaim;
    maxICAIM = -tmp_minicaim;
else
    % positive and negative
    maxICAIM = max([minICAIM, maxICAIM, threshold]);
    minICAIM = -maxICAIM;
end
