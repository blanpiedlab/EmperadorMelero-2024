%% function plotOneChannel(loc,col,cchoice,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% This function will plot one xy(z)
% Call figure before plotting it or it will plot on the last open figure
% Need to include loc table and col struct.
% Optionally call name value pairs
% color: ex 'c','r' will make points red. default color is blue
% marker size: ex 'msize',10 will make points big. default size is 1
% markershape: ex 'mshape','.' will turn markers to points. default shape is 'o'

function plotOneChannel(loc,col,varargin)
    
    % Get inputs
    p = inputParser;
    addRequired(p,'loc',@ismatrix);
    addRequired(p,'col',@isstruct);
    addParameter(p,'c','b',@ischar);
    addParameter(p,'msize',1,@isnumeric);
    addParameter(p,'mshape','o',@ischar);
    p.parse(loc,col,varargin{:});
    loc = p.Results.loc;
    col = p.Results.col;
    cchoice = p.Results.c;
    markersize = p.Results.msize;
    markershape = p.Results.mshape;

    % Plot one channel
    if isempty(col.z) == 1 % plot the 2d case
         plot(loc(:,col.x), loc(:,col.y),...
        'Color', cchoice,... %plots x and y for 647 with characteristics below in blue
        'Markersize', markersize,...
        'LineStyle', 'none',...
        'MarkerFaceColor', cchoice,...
        'Marker', markershape);
    elseif isempty(col.z) == 0 % plot the 3d case
        plot3(loc(:,col.x),loc(:,col.y),loc(:,col.z),...
        'Color', cchoice,... %plots x and y for 647 with characteristics below in blue
        'Markersize', markersize,...
        'LineStyle', 'none',...
        'MarkerFaceColor', cchoice,...
        'Marker', markershape);
    end

    axis equal
    
end
