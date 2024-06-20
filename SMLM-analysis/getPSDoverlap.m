function pctoverlap = getPSDoverlap(psdlocs,otherlocs)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Sarah Metzbower, Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
    % Inputs:
    % psdlocs - xy coordinates of a scaffold protein, eg. PSD95
    % otherlocs - xy coordinates of another other protein whose pct overlap
    % with respect to PSD95 you are interested in 
    
    %Outputs:
    % pctoverlap - in decimal between 0 and 1, the amount of PSD95's area
    % that is occupied by the other protein. Always relative to the first
    % protein's coordinates you've input. 

    p = inputParser;
    addRequired(p,'psdlocs',@isnumeric);
    addRequired(p,'otherlocs',@isnumeric);
    %addParameter(p,' ',0,@isnumeric);
    parse(p,psdlocs,otherlocs);
    
    psd = alphaShape(psdlocs(:,1), psdlocs(:,2),'HoleThreshold',100000);
    other = alphaShape(otherlocs(:,1), otherlocs(:,2),'HoleThreshold',100000);
    [~,psdboundary] = boundaryFacets(psd);
    [~,otherboundary] = boundaryFacets(other);
    % Convert to a polyshape which has more useful functions in MATLAB than
    % alphashape, including centroid and intersect
    psdshape = polyshape(psdboundary(:,1), psdboundary(:,2));
    othershape = polyshape(otherboundary(:,1), otherboundary(:,2));
    
    overlapshape = intersect(psdshape,othershape);
    
    overlaparea = area(overlapshape);
    psdarea = area(psdshape);
    pctoverlap = overlaparea/psdarea;


end
