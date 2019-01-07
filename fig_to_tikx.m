

lists = dir('fig\');

for ii = 1:length(lists)
    filename = lists(ii).name;
    if length(filename) > 4 && strcmp(filename(end-3:end), '.fig')
        if ~exist(['fig\', filename(1:end-3), 'tex'],'file')
            disp(filename)
            open(['fig\', filename]);
            matlab2tikz(['fig\', filename(1:end-3), 'tex'])
            close
        end
    end
    
end