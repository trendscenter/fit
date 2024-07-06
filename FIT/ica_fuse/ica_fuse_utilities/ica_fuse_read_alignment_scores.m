function alignments = ica_fuse_read_alignment_scores(outputDir, run_no)
%% Read alignment scores (neural net fusion)
%

%numRuns=100;
%count = 0;
%for n = 1:numRuns
fileN = fullfile(outputDir, ['alignment_run', num2str(run_no), '.txt']);
strs = textread(fileN, '%s', 'delimiter', '\n');
inds = strmatch(strs{1}, strs, 'exact');

for ii = 1:length(inds)
    % count = count + 1;
    startT = inds(ii) + 1;
    if (ii ~= length(inds))
        endT = inds(ii+1) - 1;
    else
        endT = length(strs);
    end
    
    chunks = str2num(char(strs(startT:endT)));
    
    if (ii == 1)
        alignments = zeros(size(chunks, 1), size(chunks, 2), length(inds));
    end
    
    alignments(:, :, ii) = chunks;
    
end
