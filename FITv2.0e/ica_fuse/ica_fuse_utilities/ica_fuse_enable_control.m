function ica_fuse_enable_control(tags, handles, optional)
% Enable or disable control based on the optional variable

if ~exist('optional', 'var')
    optional = 'on';
end

% % Check if the optional is neither on nor off
% if ~strcmpi(optional, 'on') &  ~strcmpi(optional, 'off')
%     optional = 'on';
% end

% disable all the controls
for ii = 1:size(tags, 1)
    currentH = findobj(handles, 'tag', deblank(tags(ii, :)));
    if ~isempty(currentH)
        set(currentH, 'enable', optional);
    end
    clear currentH;
end