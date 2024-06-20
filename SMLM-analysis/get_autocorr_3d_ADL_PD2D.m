function out = get_autocorr_3d_ADL_PD2D(A, renderpixel, rmax, psize, flag)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy, Poorna Dharmasri, Aihui Tang
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% function [out] = get_autocorr_3d_ADL(A, pixel, rmax, ps, flag)
% calculates autocorrelation function for 3 dimensional images
%
% INPUTS
% x,y,z = molecular localization coordinates in pixels
% renderpixel = pixel size to render the 3-D distribution image
% rmax = maximum r value to correlate in units of pixels
% psize = camera pixelsize conversion (nm/px)
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% normally just outputs avg_g below
% G = 3 dimensional correlation function.  x, y and z values range between
%    -rmax:rmax
% r = radius values
% avg_g = angularly averaged autocorrelation function.
% var_g = errors on angularly averaged g
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 3D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.26.10 by Sarah Veatch.
% Last updated 02.18.14 by Aihui Tang.
% Cleaned up and optimized by Aaron Levy 9/2019. Only real change was
% getting rid of the inalpha function and replacing it with
% alphaShape/inShape.

% Set flag for display
if nargin ~= 5
    flag = 0; 
end

% breakout the xyz coordinates
x = A(:,1);
y = A(:,2);
%z = A(:,3);
if isempty(x)
    out = [];
elseif ~isempty(x)
    % Set minimum xyz pos to 0, shift away from edge, and convert to
    % renderpixel space (ps/pixel)
    edge = 3;
    x = x - min(x) + edge + 1; 
    y = y - min(y) + edge + 1; 

    x = floor(x * psize / renderpixel) + 1;
    y = floor(y * psize / renderpixel) + 1;

    % Make array I1, which is a pixelated version of the synapse where the
    % value at an index is the number of pixelated locs that fall in that space
    I1 = zeros(max(x) + edge, max(y) + edge);
    for i = 1:length(x)
        I1(x(i),y(i)) = I1(x(i),y(i)) + 1;
    end

    % Make array mask that gives pixelated xyz positions of I1 that are
    % actually within the alphashape of the synapse. This first crops the mask
    % down to eliminate points between 0 and min(xyz), which is set by edge+1
    % above. This makes inShape faster.
    mask = I1;
    lim = (edge+1)*psize/renderpixel-1;     % this is the min value for xyz set above; the -1 just gives a cushion to include all points
    mask(1:lim,:) = nan;                  % points that are below the minimum xyz positions to nan
    mask(:,1:lim) = nan;
    mask(~isnan(mask)) = 0;                 % any points not nan should be 0, rather than the result of I1, for the mask
    [cx, cy] = ind2sub(size(mask), find(mask == 0)); % get the xyz indices of all points in mask that aren't nan
    marr = [cx cy]; 
    shp = alphaShape(x,y,'HoleThreshold',100000);%,1.5*psize/renderpixel);       % find the alpha shape of the synapse.
    tf = inShape(shp,marr(:,1),marr(:,2));     % logical of which rows of marr are in shp (ignoring nan).  here is where ADL replaced inalpha
    marr = marr(tf == 1,:);
    for i = 1:sum(tf)
        mask(marr(i,1),marr(i,2)) = 1;
    end
    mask(isnan(mask)) = 0;          

    % just here for seeing correct cropping of the mask before inshape
    % figure;scatter(cx,cy,'.b');xlim([0 200]);ylim([0 200]);zlim([0 200]);
    % figure;scatter(x,y,'.r');xlim([0 200]);ylim([0 200]);zlim([0 200]);

    % This is the old way for making the mask that is slower, but I'm leaving
    % here for reference.
    % Make array mask that gives the pixelated xyz positions of I1 that are
    % actually within the alpha shape of the synapse
    % mask = zeros(size(I1));
    % [cx, cy, cz] = ind2sub(size(mask), find(mask == 0));
    % marr = [cx cy cz];   % xyz are size(I1). row 1: x(1) y(1) z(1), row 2: x(2) y(2) z(2) row numel(x)+1: x(1) y(2) z(2)
    % shp = alphaShape(x,y,z,1.5*psize/renderpixel);       % find the alpha shape of the synapse.
    % tf = inShape(shp,marr(:,1),marr(:,2),marr(:,3));     % logical of which rows of marr are in shp.  here is where ADL replaced inalpha
    % marr = marr(tf == 1,:);
    % for i = 1:sum(tf)
    %     mask(marr(i,1),marr(i,2),marr(i,3)) = 1;
    % end


    % THIS BLOCK DOES THE ACTUAL AUTOCORELLATION CALCULATION
    N = sum(sum(I1.*mask));                                        % number of particles within mask
    A = sum(sum(mask));                                            % area of mask
    I1 = double(I1);                                                    % convert to double
    L1 = size(I1, 1) + rmax;                                            % size of fft2 (for zero padding)
    L2 = size(I1, 2) + rmax;                                            % size of fft2 (for zero padding)

    NP = real(fftshift(ifftn(abs(fftn(mask, size(I1)+rmax)).^2)));      % Normalization for correct boundary conditions
    tt = find(NP == 0); 
    if size(tt,1) > 0 
        NP(tt) = -1e-15;                                                % Make any 0s instead a very small value <0
    end
    G1 = A^2/N^2*real(fftshift(ifftn(abs(fftn(I1.*mask,size(I1)+rmax)).^2)))./NP; % 2D G(r) with proper normalization
    NPmask = max(0,max(max(max(NP)))/8-NP);
    G1 = G1 + 0.0001;
    G1(NPmask~=0) = 0;      % set non-relevant points equal to 0
    G = G1(floor(L1/2+1)-rmax:floor(L1/2+1)+rmax, floor(L2/2+1)-rmax:floor(L2/2+1)+rmax);  
    % G is the 3D autoC results filtered to only the valid part of G


    % Make 3d array r that is the pixel radial distance to the center xyz=0
    xvals = ones(rmax*2 + 1, rmax*2 + 1); 
    yvals = ones(rmax*2 + 1, rmax*2 + 1); 

    for i = 1:2*rmax + 1
        xvals(i,:,:) = i - rmax - 1;
        yvals(:,i,:) = i - rmax - 1;
    end
    r = sqrt(xvals.^2 + yvals.^2);

    % reshape r and G to columns, remove elements where G = 0, then
    % re-sort r and G by increasing r, and bin by radius.
    Ar = reshape(r, 1, (2*rmax+1)^2);
    Avals = reshape(G, 1, (2*rmax+1)^2);
    Ar = Ar(Avals ~= 0);
    Avals = Avals(Avals ~= 0);
    [rr,ind] = sort(Ar);                         % sort by r values
    vv = Avals(ind);                             % resort Avals to be the same order as rr
    r = 0:floor(max(rr));                        % the radii you want to extract
    [~, ~, bin] = histcounts(rr, r-.5);          % bin by radius

    % extract the average and variance of G in each bin
    avg_g = zeros(1,rmax + 1);
    var_g = zeros(1,rmax + 1);
    for j = 1:rmax + 1                                
        m = bin == j;
        n2 = sum(m);                                            % the number of pixels in that bin
        avg_g(j) = sum(m.*vv) / n2;                             % the average G values in this bin
        var_g(j) = sqrt(sum(m.*(vv - avg_g(j)).^2)) / n2;       % the variance of the mean
    end

    G(rmax+1, rmax+1, rmax+1) = 0;

    if flag
        r = [0:rmax].*renderpixel;
        figure('Color','white')
        errorbar(r(2:length(r)), avg_g(2:length(r)), var_g(2:length(r)));
        axis tight
    end

    out = avg_g; %[r, avg_g, var_g];
end
end % function
