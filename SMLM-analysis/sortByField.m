%% out = sortByField(in,sortcol)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
% sorts an input loc table loaded from Omniloader (in) by the
% column sortcol into a cell array
function out = sortByField(in,sortcol)

    sorted = sortrows(in,sortcol);
    tabled = tabulate(sorted(:,sortcol)); % get frequency table of how many times each frame appears
    out = mat2cell(sorted,tabled(:,2));
    
end
