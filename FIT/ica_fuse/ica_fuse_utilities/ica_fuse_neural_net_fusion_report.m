function ica_fuse_neural_net_fusion_report(param_file, results)
%% Report generator for neural net fusion
%

sesDir = fileparts(param_file);

if (isempty(sesDir))
    sesDir = pwd;
end

load (param_file);

try
    formatName = results.formatName;
catch
end

if (~exist('formatName', 'var'))
    formatName = 'html';
end

if (isnumeric(formatName))
    if (formatName == 1)
        formatName = 'html';
    else
        formatName = 'pdf';
    end
end

results.formatName = formatName;


resultsFile = 'ica_fuse_stats_neural_net_fusion';
outDir = fullfile(sesDir, [neuralNetFusionInfo.prefix, '_neural_net_results']);
opts.codeToEvaluate = 'ica_fuse_stats_neural_net_fusion(param_file, results);';


results.outputDir = outDir;
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opts.useNewFigure = false;
end


disp('Generating reults summary. Please wait ....');
drawnow;

assignin('base', 'param_file', param_file);
assignin('base', 'results', results);


publish(resultsFile, opts);


if (strcmpi(opts.format, 'html'))
    ica_fuse_openHTMLHelpFile(fullfile(outDir, [resultsFile, '.html']));
else
    open(fullfile(outDir, [resultsFile, '.pdf']));
end

disp('Done');


close all;