%% function out = get_crosscorr_2d(A, B, psize, renderpixel, rmax, varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy, Aihui Tang, Poorna Dharmasri
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% calculates crosscorrelation function for 2 dimensional images
% PLEASE NOTE that this version does NOT do any shifting of the data
% relative to the peak, as in the 3d code (ie xyzshift). For 2d data, we
% don't feel it is appropriate to shift the data based on the smoothed
% cross correlation. If that were to be implemented, you would need to work
% out better smoothing/density smashing here so the shift is not chasing noise as it would
% be now, and will only fit the overall shape - ADL and PD 5/17/23
%
% REQUIRED INPUTS
% A = x and y points only (ie m x 2 array) of channel 2
% B = same as A for channel 2
% psize = pixelsize of the image in nm/px
% rmax = maximum r value to correlate in units of pixels, usually 50
% (deprecated) maxdistance = max distance over which to search for peaks, often 180 but
% sometimes less
% 
% OPTIONAL INPUTS
% figflag = defaults 0, set to 1 to display figure of xcorr
% smooth = default 1, set to >1 (ie 3) to smooth the cross correlation and
% suppress noise. be careful though if you go too high it can suppress the
% xc entirely
%
% OUTPUTS
% out = angularly averaged autocorrelation function. NOTE that usually we
% suppress the first point.
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 3D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.26.10 by Sarah Veatch.
% Last updated 02.18.14 by Aihui Tang.
% Last updated 5/1/23 by ADL for cleanup, alphashape update, function
% Last updated 5/17/23 ADL/PD to get rid of the shift for peak and
% maxdistance as input
% inputs, turn off smoothing default

function out = get_crosscorr_2d_ADL(A, B, psize, renderpixel, rmax, varargin)

    p = inputParser;
    addRequired(p,'A',@isnumeric)
    addRequired(p,'B',@isnumeric)
    addRequired(p,'psize',@isnumeric)
    addRequired(p,'renderpixel',@isnumeric)
    addRequired(p,'rmax',@isnumeric)
    %addRequired(p,'maxdistance',@isnumeric)
    addParameter(p,'figflag',0,@isnumeric)
    addParameter(p,'smooth',1,@isnumeric)

    parse(p,A,B,psize,renderpixel,rmax,varargin{:})
    A = p.Results.A;
    B = p.Results.B;
    psize = p.Results.psize;
    renderpixel = p.Results.renderpixel;
    rmax = p.Results.rmax;
    %maxdistance = p.Results.maxdistance;
    figflag = p.Results.figflag;
    smoothval = p.Results.smooth;

    % do some data checks
    if size(A,2)~=2 || size(B,2)~=2
        error('Dimensions of A or B is wrong; only input xy data')
    end

    %%%% Make a pixelized render vased on renderpixel
    % Get x and y coordinates
    x1 = A(:,1);
    y1 = A(:,2); 
    x2 = B(:,1); 
    y2 = B(:,2); 

    % Zero and shift over by edge + 1
    edge=3;
    x0 = min(min([x1;x2])); 
    y0 = min(min([y1;y2])); 
    x1 = x1-x0+edge+1; 
    y1 = y1-y0+edge+1; 
    x2 = x2-x0+edge+1; 
    y2 = y2-y0+edge+1; 

    % Convert to renderpixels
    x1 = floor(x1*psize/renderpixel); 
    x2 = floor(x2*psize/renderpixel);
    y1 = floor(y1*psize/renderpixel); 
    y2 = floor(y2*psize/renderpixel);

    % Make I1 and I2, which are pixel versions of synapses
    x0 = max(max([x1;x2]))+edge*psize/renderpixel; 
    y0 = max(max([y1;y2]))+edge*psize/renderpixel; 

    I1 = zeros(x0,y0);
    for i = 1:length(x1)
        I1(x1(i),y1(i)) = I1(x1(i),y1(i)) + 1;
    end
    I2 = zeros(x0,y0);
    for i = 1:length(x2)
        I2(x2(i),y2(i)) = I1(x2(i),y2(i)) + 1;
    end

    warning('off','MATLAB:alphaShape:DupPointsBasicWarnId');
    % generate masks for both channels
    mask1 = I1;
    lim = min([x1;y1])-1;                 % this is the min value for xyz set above with some safety
    mask1(1:lim,:) = nan;                  % set points that are below the minimum xyz positions to nan
    mask1(:,1:lim) = nan;
    mask1(~isnan(mask1)) = 0;  
    [cx, cy] = ind2sub(size(mask1), find(mask1 == 0));
    marr = [cx cy]; 
    shp = alphaShape(x1,y1,'HoleThreshold',100000);            % create alpha shape for synapse %formerly 1.5*psize/renderpixel but PAD changed on 01/30/23 to 150/renderpixel to preserve proper radius per Aihui OG analysis
    tf = inShape(shp,marr(:,1),marr(:,2));             % find pixelated positions that are in the alpha shape
    marr = marr(tf == 1,:);
    for i=1:sum(tf)
        mask1(marr(i,1),marr(i,2)) = 1;
    end
    mask1(isnan(mask1)) = 0;

    mask2 = I2;
    lim2 = min([x2;y2])-1;                 % this is the min value for xyz set above with some safety
    mask2(1:lim2,:) = nan;                  % set points that are below the minimum xyz positions to nan
    mask2(:,1:lim2) = nan;
    mask2(~isnan(mask2)) = 0;  
    [cx2, cy2] = ind2sub(size(mask2), find(mask2 == 0));
    marr2 = [cx2 cy2]; 
    shp2 = alphaShape(x2,y2,'HoleThreshold',100000);            % create alpha shape for synapse %formerly 1.5*psize/renderpixel but PAD changed on 01/30/23 to 150/renderpixel to preserve proper radius per Aihui OG analysis
    tf2 = inShape(shp2,marr2(:,1),marr2(:,2));             % find pixelated positions that are in the alpha shape
    marr2 = marr2(tf2 == 1,:);
    for i=1:sum(tf2)
        mask2(marr2(i,1),marr2(i,2)) = 1;
    end
    mask2(isnan(mask2)) = 0;
    warning('on','MATLAB:alphaShape:DupPointsBasicWarnId');

    %Do the xcorr

    N1 = sum(sum(sum(I1.*mask1)));  % number of particles within channel 1
    N2 = sum(sum(sum(I2.*mask2)));  % number of particles within channel 2
    A1 = sum(sum(sum(mask1)));      % area of mask1
    A2 = sum(sum(sum(mask2)));      % area of mask2
    I1 = double(I1);  
    I2 = double(I2);       % convert to double

    NP = real(fftshift(ifftn(fftn(mask1,size(I1)+rmax).*conj(fftn(mask2,size(I1)+rmax)))));
    tt = find(NP==0); 
    if size(tt,1)>0
        NP(tt)=-1e-15;
    end

    G1 = A1*A2/N1/N2*real(fftshift(ifftn(fftn(I1.*mask1,size(I1)+rmax).*conj(fftn(I2.*mask2,size(I1)+rmax)))))./NP; % 2D G(r) with proper normalization
    NPmask = max(0,max(max(NP))/4-NP);
    G1(NPmask ~= 0)=0; 

    % smooth the data (if smoothval = 1 there is no smoothing)
    k = ones(smoothval,smoothval)./smoothval^2;
    G1 = conv2(abs(G1),k,'same');

    % look for peaks within md of the center of the xcorr
%     G = G1*0+1;
%     [x,y] = ind2sub(size(G), find(G));
     L = size(G1);
     L = floor(L/2+1);
%     d = pdist2([x, y],L);       % measure distance of each xy position to center
%     md = maxdistance / renderpixel;   % allowable distance to look for peak
%     G(d>md) = 0;                % set anything outside distance to 0
%     G = G.*G1;                  % mask G1 by resulting window
%     Gm = max(max(max(G)));      % get the max value
%     [mx, my] = ind2sub(size(G), find(G==Gm));         % find the position of the peak of correlation
%     G = G1(mx-rmax:mx+rmax, my-rmax:my+rmax);         % crop G1 to within rmax of the peak
    G = G1(L(1)-rmax:L(1)+rmax, L(2)-rmax:L(2)+rmax);         % crop G1 to within rmax of the peak

    % if figflag
    %     figure('Color', 'white'); imshow(abs(G1))
    %     figure('Color', 'white'); imshow(abs(G))
    % end

    xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    % map to x positions with center x=0
    yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    % map to y positions with center y=0
    zvals = G;
    [~,r,v] = cart2pol(xvals,yvals, zvals);     % convert x, y to polar coordinates

    Ar = reshape(r,1, (2*rmax+1)^2);
    Avals = reshape(v, 1, (2*rmax+1)^2);
    temp = find(Avals ~=0);
    Ar = Ar(temp); 
    Avals = Avals(temp);
    [rr,ind] = sort(Ar);                         % sort by r values
    vv = Avals(ind);                             % reindex g
    r = 0:floor(max(rr));                        % the radii you want to extract
    [~,~, bin] = histcounts(rr, r-.5);           % bin by radius

    % this is the old strategy for posterity; below is the same with
    % preallocated arrays
    % for j = 1:rmax+1                             % now get averages
    %     m = bin==j;
    %     n2 = sum(m);                             % the number of pixels in that bin
    %     if n2==0
    %         vals(j)=0; 
    %         er(j)=0;                             % if no bins, no data
    %     else
    %         g(j) = sum(m.*vv)/n2;                % the average G values in this bin
    %         dg(j) = sqrt(sum(m.*(vv-g(j)).^2))/n2; % the variance of the mean
    %     end
    % end

    % zero values will leave nan instead of zero when there are no values
    g = nan(1,rmax + 1);
    dg = nan(1,rmax + 1);
    vals = nan(1,rmax + 1);
    er = nan(1,rmax + 1);
    for j = 1:rmax + 1                                
        m = bin == j;
        n2 = sum(m);                                            % the number of pixels in that bin
        if n2==0
            vals(j) = 0;
            er(j) = 0;
        else
            g(j) = sum(m.*vv) / n2;                             % the average G values in this bin
            dg(j) = sqrt(sum(m.*(vv - g(j)).^2)) / n2;       % the variance of the mean
        end
    end

    %G(rmax+1, rmax+1) = 0; % ADL here- i dont know why this is here but it doesn't do anything.

    if figflag
        r = 0:rmax;
        figure('Color', 'white');
        errorbar(r(1:length(r)), g(1:length(r)), dg(1:length(r)));
        axis tight
    end

    out = g; %[r; g; dg];
    
end

