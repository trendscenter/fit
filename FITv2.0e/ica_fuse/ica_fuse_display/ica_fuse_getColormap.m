function cm = ica_fuse_getColormap(imageValues, numOfComp, displayType)
% Get the colormap based on the image values
% Colors are obtained from colors.mat file

ica_fuse_defaults;
global COLORMAPFILE;

if ~exist('displayType', 'var')
    displayType = 'other';
end

load(COLORMAPFILE);

%--set up colormap
colorLength = 64;
for i=1:numOfComp
    %get colormaps
    if(i ==1)
        if strcmp(lower(displayType), 'composite')
            cm = redV;
        else
            if(imageValues == 2 | imageValues ==3)
                cm = hot2;
            elseif(imageValues == 4) % colormap for negative values
                cm = cold;
                % flip the color range for negative image values
                [nrows, ncols] = size(cold);
                for ii = 1: nrows
                    cm(ii, 1:ncols) = cold(nrows - ii + 1, 1:ncols);
                end
            else
                cm = coldhot2; %coldhot;
            end
        end
    elseif(i==2)
        if(imageValues == 2 | imageValues ==3)
            cm = blueV; %green;
        else
            cm = blueV; %green;
        end
    elseif(i==3)
        if(imageValues == 2 | imageValues ==3)
            cm = greenV; %red;
        else
            cm = greenV; %red;
        end
    elseif(i==4)
        if(imageValues == 2 | imageValues ==3)
            cm = pinkV; %purple;
        else
            cm = pinkV; %purple;
        end
    elseif(i==5)
        if(imageValues == 2 | imageValues ==3)
            cm = yellowV; %blue;
        else
            cm = yellowV; %blue;
        end
    else
        if(imageValues == 2 | imageValues ==3)
            cm = yellowV; %orange;
        else
            cm = yellowV; %orange;
        end
    end

    colorbarSkip = size(cm,1)/colorLength;
    if(i==1)
        tempCM = [cm(1:colorbarSkip:end,:)];
    else
        tempCM = [tempCM;cm(1:colorbarSkip:end,:)];
    end
end
cm = tempCM;

gm = linspace(0, 1, 64); 
if size(gm, 1) == 1
    gm = gm';
end
structCM = [gm gm gm];
colorbarSkip = size(structCM,1)/colorLength;
cm = [cm; structCM(1:colorbarSkip:end, :)];

colorbarSkip = ceil(size(cm,1)/256);
cm = cm(1:colorbarSkip:end, :);