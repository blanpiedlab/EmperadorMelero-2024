%% Centrality - Does a test position within an elliptic synapse occur closer to the center or edge of the ellipse?
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Poorna Dharmasri
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Written by PAD on 05/08/2023
% REQUIRES
%  REQURIED INPUTS-
%
%   clusterlocs: the xy coordinates of the localizations that define the
%   region describing the synapse (e.g. PSD-95 locs).
%
%   testloc: an nx2 array containing the xy coordinate(s) of the loc
%   whose position within the synapse you want to test.
%   NOTE: Any number of locs can be passed at once
%   (e.g. just one loc if you're looking at a nc center or an array of
%   NC centers or every loc in the synapse, it can handle it all).
%
%  OPTIONAL INPUTS
%   plotflag: default 0, use as Name Value with value 1 to plot a figure
%   that shows the full synaptic ellipse, the position of the test loc,
%   and, depending on the analysis, the relevant graph to show the
%   component parts of the outcome of that analysis
%
%   legacy: default 0, use as Name Value pair with value 1 to use the
%   legacy version of this code. This version uses a different method to
%   calculate the centrality which involves the equation of the ellipse
%   modified to include a rotation around the z axis by an angle theta.
%   This yields the same result as the default linear algebra approach, but
%   is computationally slower in matlab.
%
%   unitcircle: default 0, use as Name Value with value 1 to calculate a
%   radial position of the testloc as a fractional area of a unit circle
%   (as in Li et al., 2021 from Watanabe lab).
%
% OUTPUTS-
%   centrality: Represents the fractional area occupied by an ellipse that
%   contains the testloc on its perimeter, relative to the ellipse that 
%   describes the synapse. A value of 0.25 is perfectly centered - less 
%   than that is closer to center, greater than that is closer to edge
%
%   varargout{1}: The sqrt of the centrality output, this essentially is a
%   normalized radial position of the testloc that takes into account the
%   potential difference in length of the semi-minor and semi-major axes.
%   It essentially converts the readout of this analysis into a linear space, ie.
%   values of 0-0.5 are closer to center of synapse, 0.5-1 are closer to
%   edge.
%
%   varargout{2}: The unit cirlce area - based off of Li et al., 2021, assuming
%   a unit circle using a vector from center of synapse to test loc to
%   intersection with synaptic border beyond, this is the squared
%   fractional distance along that vector that the testloc is found on.
%
%  WHAT THE FUNCTION DOES:
%  This function treats all points in the synapse as if they belong to a
%  coordinate system that is described by the semimajor and semiminor axes
%  of the ellipse. That is, the xy coordinate of the testloc is not really
%  described by xy, but rather by x'y' relative to the axes of the ellipse.
%  Thus, to use the ellipse equation to understand the fractional area, we
%  must first shift the testloc into xy space. To do that, we first 
%  calculate the basis vectors that describe the coordinate system upon 
%  which the ellipse fit around the synapse lies, ie. accounts for rotation
%  around the z axis. We scale the vectors from the center of the ellipse 
%  to each vertex encompassing the testloc such that their magnitude is 
%  equal to 1. These are concatenated into a transformation matrix, as they
%  describe the transformation that is applied to the standard basis vectors to
%  arrive at the new coordinate system. The inverse of this transformation
%  matrix is then calculated and then applied to the original coordinates 
%  of the testlocs. Thus, the original testloc, which belonged is x'y' space
%  is shifted into xy space, allowing us to use the equation of an ellipse 
%  to then calculate the fractional area of an ellipse whose perimeter
%  contains the transformed xy point. 


function [centrality,varargout] = Centrality(clusterlocs,testloc,varargin)


    p = inputParser;
    addRequired(p,'clusterlocs');
    addRequired(p,'testloc');
    addParameter(p,'plotflag',0,@isnumeric);
    addParameter(p,'legacy',0,@isnumeric);
    addParameter(p,'unitcircle',0,@isnumeric);
    parse(p,clusterlocs,testloc,varargin{:});

    %Lines 46 through 152 are from Nima Moshtagh (2023). Minimum Volume Enclosing Ellipsoid (https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid), MATLAB Central File Exchange. Retrieved May 4, 2023.
    %with minor modifications by PAD
    P = p.Results.clusterlocs';
    tolerance = 0.01;

    % ---------------------------------
    % data points
    % -----------------------------------
    [d, N] = size(P);
    Q = zeros(d+1,N);
    Q(1:d,:) = P(1:d,1:N);
    Q(d+1,:) = ones(1,N);
    % initializations
    % -----------------------------------
    count = 1;
    err = 1;
    u = (1/N) * ones(N,1);          % 1st iteration
    % Khachiyan Algorithm
    % -----------------------------------
    while err > tolerance
        X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix
        M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
        [maximum, j] = max(M);
        step_size = (maximum - d -1)/((d+1)*(maximum-1));
        new_u = (1 - step_size)*u ;
        new_u(j) = new_u(j) + step_size;
        count = count + 1;
        err = norm(new_u - u);
        u = new_u;
    end
    %%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
    % Finds the ellipse equation in the 'center form':
    % (x-c)' * A * (x-c) = 1
    % It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
    % of the ellipse.
    U = diag(u);
    % the A matrix for the ellipse
    % --------------------------------------------
    A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );
    % center of the ellipse
    % --------------------------------------------
    ellipsecenter = P * u;

    C = ellipsecenter;
    N = .5; % Default value for grid
    % See if the user wants a different value for N.
    %----------------------------------------------
    %     if nargin > 2
    %         N = varargin{1};
    %     end
    % check the dimension of the inputs: 2D or 3D
    %--------------------------------------------
    if length(C) == 3
        Type = '3D';
    elseif length(C) == 2
        Type = '2D';
    else
        disp('Cannot plot an ellipse with more than 3 dimensions!!');
        return
    end
    % "singular value decomposition" to extract the orientation and the
    % axes of the ellipsoid
    [U, D, V] = svd(A);
    if strcmp(Type, '2D')
        % get the major and minor axes
        %------------------------------------
        a = 1/sqrt(D(1,1));
        b = 1/sqrt(D(2,2));
        axis1radius = a;
        axis2radius = b;
        theta = [0:N:359.5].*(pi/180);%[0:1/N:2*pi+1/N]; % switched to 0:N:359.5 after some errors in Kaeser project where 
        % first point was duplicated in problematic way. Seems to solve the
        % problem here.
        % Parametric equation of the ellipse
        %----------------------------------------
        state(1,:) = a*cos(theta);
        state(2,:) = b*sin(theta);
        axis1pointinds = [find(state(1,:)==max(state(1,:))) find(state(1,:)==min(state(1,:)))];
        axis2pointinds = [find(state(2,:)==max(state(2,:))) find(state(2,:)==min(state(2,:)))];
        % Coordinate transform
        %----------------------------------------
        X = V * state;
        X(1,:) = X(1,:) + C(1);
        X(2,:) = X(2,:) + C(2);
        axis1pointsall = unique(X(:,axis1pointinds)','rows');
        axis2pointsall = unique(X(:,axis2pointinds)','rows');
        ellipsecenter = ellipsecenter';
        ellipseshape = polyshape(X');
    elseif strcmp(Type,'3D')
        % generate the ellipsoid at (0,0,0)
        %----------------------------------
        a = 1/sqrt(D(1,1));
        b = 1/sqrt(D(2,2));
        ellipsecenter = 1/sqrt(D(3,3));
        [X,Y,Z] = ellipsoid(0,0,0,a,b,ellipsecenter,N);

        %  rotate and center the ellipsoid to the actual center point
        %------------------------------------------------------------
        XX = zeros(N+1,N+1);
        YY = zeros(N+1,N+1);
        ZZ = zeros(N+1,N+1);
        for k = 1:length(X)
            for j = 1:length(X)
                point = [X(k,j) Y(k,j) Z(k,j)]';
                P = V * point;
                XX(k,j) = P(1)+C(1);
                YY(k,j) = P(2)+C(2);
                ZZ(k,j) = P(3)+C(3);
            end
        end
    end

    if p.Results.legacy == 0
        findmajoraxispos = axis2pointsall-ellipsecenter;
        tomajorbasis = axis2pointsall(findmajoraxispos(:,1)>=0,:);

        findminoraxispos = axis1pointsall-ellipsecenter;
        tominorbasis = axis1pointsall(findminoraxispos(:,2)>=0,:);

        majorbasis = (tomajorbasis'-ellipsecenter');
        majorbasis = majorbasis/norm(majorbasis);
        minorbasis = tominorbasis'-ellipsecenter';
        minorbasis = minorbasis/norm(minorbasis);

        transform_matrix = [majorbasis minorbasis];
        OGbasis_testlocvector = p.Results.testloc'-ellipsecenter';
        inverse_tm = inv(transform_matrix);
        new_xy = inverse_tm*OGbasis_testlocvector;
        tocentrality = ((new_xy(1,:).^2)./(axis2radius^2)) + ((new_xy(2,:).^2)./(axis1radius^2));
        centrality = tocentrality';
        varargout{1} = sqrt(centrality);

    elseif p.Results.legacy == 1
        dist_to_a1 = pdist2(axis1pointsall,p.Results.testloc);
        a1 = cell(size(dist_to_a1,2),1);
        for i = 1:size(dist_to_a1,2)
            a1{i,1} = [ellipsecenter; axis1pointsall(dist_to_a1(:,i)==min(dist_to_a1(:,i)),:)];
        end

        dist_to_a2 = pdist2(axis2pointsall,p.Results.testloc);
        a2 = cell(size(dist_to_a2,2),1);
        for i = 1:size(dist_to_a2,2)
            a2{i,1} = [ellipsecenter; axis2pointsall(dist_to_a2(:,i)==min(dist_to_a2(:,i)),:)];
        end

        x1 = p.Results.testloc(:,1);
        y1 = p.Results.testloc(:,2);
        h = ellipsecenter(1);
        k = ellipsecenter(2);
        majoraxesvectors = axis2pointsall-[h k];
        a = norm(majoraxesvectors(1,:));
        minoraxesvectors = axis1pointsall-[h k];
        b = norm(minoraxesvectors(1,:));
        xaxisvector = [h+norm(majoraxesvectors(1,:)) k]- [h k];
        thismajoraxisvector = majoraxesvectors(majoraxesvectors(:,2)>0,:);
        theta = acos(dot(thismajoraxisvector,xaxisvector)/(norm(thismajoraxisvector)*norm(xaxisvector)));

        term1 = (((x1-h)*cos(theta)) + ((y1-k)*sin(theta))).^2;
        term2 = (((((x1-h)*sin(theta)) - ((y1-k)*cos(theta))).^2)./((b^2)/(a^2)));
        anew = sqrt(term1+term2);
        bnew = anew*(b/a);

        centrality = ((anew./a).*(bnew./b));
        varargout{1} = sqrt(centrality);
    end

    if p.Results.unitcircle
        vectortonccenter = p.Results.testloc-ellipsecenter;
        normalizedvector = vectortonccenter./vecnorm(vectortonccenter,2,2);
        linepastshape = ellipsecenter+(10*normalizedvector);
        extendedncvector = cell(size(linepastshape,1),1);
        for i = 1:size(linepastshape,1)
            extendedncvector{i,1} = [ellipsecenter;linepastshape(i,:)];
        end
        intersectionwithellipse = cellfun(@(x) intersect(ellipseshape,x),extendedncvector,'uni',0);
        dist_to_testloc = pdist2(ellipsecenter,p.Results.testloc)';
        length_of_intersection = cellfun(@(x) pdist2(x(1,:),x(2,:)),intersectionwithellipse,'uni',0);
        unitcirclearea = (dist_to_testloc./cell2mat(length_of_intersection)).^2;

        if p.Results.plotflag
            th = 0:pi/50:2*pi;
            unitcircle = cell(size(p.Results.testloc,1),1);
            for thisloc = 1:size(p.Results.testloc,1)
                xunit = pdist2(ellipsecenter,p.Results.testloc(thisloc,:)) * cos(th) + ellipsecenter(1);
                yunit = pdist2(ellipsecenter,p.Results.testloc(thisloc,:)) * sin(th) + ellipsecenter(2);
                unitcircle{thisloc,1} = polyshape(xunit,yunit);
            end
        end
        varargout{2} = unitcirclearea;
    end

    if p.Results.plotflag && p.Results.legacy == 0
        for fignum = 1:size(p.Results.testloc)
            figure
            hold on;
            axis equal
            plot(ellipseshape)
            if p.Results.unitcircle
                plot(unitcircle{fignum})
            end
            scatter(ellipsecenter(1),ellipsecenter(2),'ok');
            plot([ellipsecenter(1);majorbasis(1)+ellipsecenter(1)], [ellipsecenter(2);majorbasis(2)+ellipsecenter(2)],'r');
            plot([ellipsecenter(1);minorbasis(1)+ellipsecenter(1)], [ellipsecenter(2);minorbasis(2)+ellipsecenter(2)],'r');
            plot([ellipsecenter(1);1+ellipsecenter(1)], [ellipsecenter(2);0+ellipsecenter(2)],'b');
            plot([ellipsecenter(1);0+ellipsecenter(1)], [ellipsecenter(2);1+ellipsecenter(2)],'b');
            scatter(p.Results.clusterlocs(:,1), p.Results.clusterlocs(:,2),'.k');
            scatter(axis1pointsall(:,1), axis1pointsall(:,2),'og');
            scatter(axis2pointsall(:,1), axis2pointsall(:,2),'og');
            scatter(p.Results.testloc(fignum,1),p.Results.testloc(fignum,2),'xr','LineWidth',3);
            scatter(new_xy(1,fignum)+ellipsecenter(1),new_xy(2,fignum)+ellipsecenter(2),'xb','LineWidth',3);
        end
    elseif p.Results.plotflag && p.Results.legacy == 1

        newellipseshapes = cell(size(p.Results.testloc,1),1);
        for thisloc = 1:size(p.Results.testloc,1)
            xplot = [];
            yplot = [];
            for rad = 0:1/20:2*pi+1/20
                xplot = [xplot;(anew(thisloc)*cos(rad)*cos(theta))-(bnew(thisloc)*sin(rad)*sin(theta))+h];
                yplot = [yplot;(anew(thisloc)*cos(rad)*sin(theta))+(bnew(thisloc)*sin(rad)*cos(theta))+k];
            end
            newellipseshapes{thisloc,1} = [xplot yplot];
        end

        vectortoa1 = cell2mat(cellfun(@(x) x(2,:)-x(1,:),a1,'uni',0));
        a1vectormag = vectortoa1./vecnorm(vectortoa1,2,2);
        bnewpoints = ellipsecenter + (bnew.*a1vectormag);

        vectortoa2 = cell2mat(cellfun(@(x) x(2,:)-x(1,:),a2,'uni',0));
        a2vectormag = vectortoa2./vecnorm(vectortoa2,2,2);
        anewpoints = ellipsecenter + (anew.*a2vectormag);

        for fignum = 1:size(p.Results.testloc)
            figure
            hold on;
            axis equal
            plot(ellipseshape)
            if p.Results.unitcircle
                plot(unitcircle{fignum})
            end
            scatter(ellipsecenter(1),ellipsecenter(2),'ok');
            scatter(a1{fignum}(2,1),a1{fignum}(2,2),'xg');
            scatter(a2{fignum}(2,1),a2{fignum}(2,2),'xg');
            plot(a1{fignum}(:,1),a1{fignum}(:,2),'k');
            plot(a2{fignum}(:,1),a2{fignum}(:,2),'k');
            plot([h;anewpoints(fignum,1)],[k;anewpoints(fignum,2)],'r');
            plot([h;bnewpoints(fignum,1)],[k;bnewpoints(fignum,2)],'r');
            scatter(p.Results.clusterlocs(:,1), p.Results.clusterlocs(:,2),'.k');
            scatter(p.Results.testloc(fignum,1),p.Results.testloc(fignum,2),'xb','LineWidth',3);
            plot(newellipseshapes{fignum}(:,1),newellipseshapes{fignum}(:,2),'r')
            if ~p.Results.unitcircle
                title(['Center to Edge measure result: ',num2str(centrality(fignum))]);
            elseif p.Results.unitcircle
                title(['Center to Edge measure result: ',num2str(centrality(fignum)), ' And unit circle result: ',num2str(unitcirclearea(fignum))]);
            end
        end
    end

end



