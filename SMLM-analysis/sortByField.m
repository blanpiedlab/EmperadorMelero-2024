%% out = sortByField(in,sortcol)
% sorts an input loc table loaded from Omniloader (in) by the
% column sortcol into a cell array
function out = sortByField(in,sortcol)

    sorted = sortrows(in,sortcol);
    tabled = tabulate(sorted(:,sortcol)); % get frequency table of how many times each frame appears
    out = mat2cell(sorted,tabled(:,2));
    
end