function [ratio, longaxis, shortaxis] = getPSDaxes_alphaRad_ADL(xycoordinates,varargin)
    %This function pulls out the long and short axes of an irregular shape,
    %such as the border of a synapse. It utilizes x and y coordinates (and
    %thus is suited only for 2D synapses as of this writing) to plot the
    %synapse and get its border, from which a centroid is determined. I
    %then create radial lines out from the center that are separated by 0.5
    %degree angles to appropriately sample the border of the psd. I then
    % use this to isolate lines that intersect the borders of the psd and
    % determine which of these lines are the longest and shortest. 
    %
    %
    % The inputs of this function are: 
    % xycoordinates: IN PIXELS, the xycoordinates of the psd you want
    % analyzed.
    % OPTIONAL: plotaxes. Default is 0, but when toggled to 1 will plot the
    % shape of the PSD, the long axis, the short axis, and the centroid.
    % alphaRadius: default is 0 and alpha radius will be minimal. when
    % toggled to any other value will calculate alpha shape with that alpha
    % radius, suggested is 2.
    %
    % The outputs of this function are: 
    % ratio = the length of the long axis divided by the length of the
    % short axis
    % longaxis = the xy coordinates of the long axis line
    % shortaxis = the xy coordinates of the short axis line
     
    %Written by PD on 03.31.22 for NMDAR project
    
    %Parse inputs
    p = inputParser;
    addRequired(p,'xycoordinates',@isnumeric);
    addParameter(p,'plotaxes',0,@isnumeric);
    addParameter(p,'alphaRadius',0,@isnumeric)
    parse(p,xycoordinates,varargin{:});


    % Use an alphashape of the PSD, drawn with some buffer space (the
    % factor 2 at the end) to get the border of the PSD
    if p.Results.alphaRadius == 0
        testshape = alphaShape(xycoordinates(:,1), xycoordinates(:,2),'HoleThreshold',100000);
    elseif p.Results.alphaRadius ~= 0
        testshape = alphaShape(xycoordinates(:,1), xycoordinates(:,2),p.Results.alphaRadius);
    end
    % Get the vertices of the synaptic border
    [~,psdboundary] = boundaryFacets(testshape);
    % Convert to a polyshape which has more useful functions in MATLAB than
    % alphashape, including centroid and intersect
    psdshape = polyshape(psdboundary(:,1), psdboundary(:,2));
    % Get the centroid
    [cx,cy] = centroid(psdshape);
    % Get the distances of each boundary vertex from the centroid. We will
    % use this to get the maximum vertex point (max(distovert)) and then
    % extend that distance by 1 px
    disttovert = pdist2([cx,cy],psdboundary);
    % Here we draw our radial lines. I may come back later and make the
    % sampling a variable input, but for now it draws lines from the center
    % that are shifted by 0.5 degrees around a circle. The x and y
    % coordinates of the distal point of the line segment are stored in
    % radlinecell
    radlinecell = {};
    anglerep = 0:.5:359.5;
    for i = 1:size(anglerep,2)
        angle = anglerep(i)*(pi/180);
        x=cx+((max(disttovert)+1)*cos(angle));
        y=cy+((max(disttovert)+1)*sin(angle));
        radlinecell{i} = [x,y];
    end

    %Here we get the lengths of all the segments that are bounded by the
    %psd shape. 
    inseglengths = NaN(360,1);
    insegs = {};
    % We only need to go halfway around the circle, since we use the
    % opposite point on the circle to complete the line segment we are
    % testing. Hence the full range of i in this loop needs to be half of
    % the max size of radlinecell. 
    for i = 1:360
        %Generate a line segment that goes through the centroid and uses
        %points i and i+360 (the opposing point on the circle)
        lineseg = [radlinecell{i}; cx cy; radlinecell{i+360}];
        %Use intersect function to ouput the coordinates along that
        %linesegment that are only in the psd. 
        [in,~] = intersect(psdshape,lineseg);
        %Store the linesegment inside the psd for later
        insegs{i} = in;
        %Get the length of that line segment by skipping the centroid point
        inlength = pdist2(in(1,:),in(3,:));
        %Store the length for each line
        inseglengths(i) = inlength;
    end

    %Get the long axis length and short axis length and use it to calculate
    %the ratio.
    lalength = max(inseglengths);
    salength = min(inseglengths);
    ratio = lalength/salength;
    % Pull the coordinates of the long and short axis
    longaxis = insegs{inseglengths==lalength};
    shortaxis = insegs{inseglengths==salength};

    %Plot the axses, centroid, and shape if needed.
    if p.Results.plotaxes == 1
        figure;
        axis equal;
        hold on;
        plot(psdshape);
        plot(longaxis(:,1), longaxis(:,2))
        plot(shortaxis(:,1), shortaxis(:,2))
        scatter(cx,cy, 'Marker', 'x', 'MarkerEdgeColor', 'r')
        title(['Ratio is ' num2str(ratio)])
    end

end

