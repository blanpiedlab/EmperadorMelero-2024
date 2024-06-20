function xyzshift = PAINTshift(A, B, renderpixel, psize, rmax, xyzshift, distance)
% THIS VERSION of the cross correlation shift is specified for re-aligning
% PAINT acquisitions with gold nanoparticles in the field after DUV
% correction. In theory, the gold will dominate the cross correlation and
% allow for any correction of a linear shift between fields leftover after
% DUV correction. It's adapted from our usual 2D xcorr code to be faster
% and get rid of some functionality we don't need here.
% 
%       For two clusters supposed to be overlaped, set xyzshift=[0,0,0];
%       For two clusters not overlaped, set xyzshift=[] if the best
%       overlaping vector is unknown. Use distance=[min, max] to limit the
%       |xyzshift|. The function will find the vector xyzshift to get the
%       best overlap between min(A,Amax/mv) and min(B,Bmax/mv). 
%
% INPUTS
% A = molecular localization coordinates (x,y,z) in pixels of channel 1
% B = molecular localization coordinates (x,y,z) in pixels of channel 2
% renderpixel = pixel size to render the 3-D distribution image, in nm; e.g. 5
% psize = camera pixel conversion (nm/px)
% rmax = maximum r value to correlate in units of pixels. default is 100;
% xyzshift = xyz vector to overlap the two clusters.
% distance = [min, max], distance variation allowed for xyzshift. input in pixels.

%
% OUTPUTS
% [xyzshift sqrt(sum(xyzshift.^2)) 0 g];
% xyzshift is the [x y z] shift to max overlap of A and B
% sqrt is essentially the MSD of xyzshift
% g is the average crosscorrelation across bins
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 3D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Original 2d xcorr code Updated 01.26.10 by Sarah Veatch.
% Modified 02.18.14 by Aihui Tang.
% Last updated 03.02.17 by Aihui Tang.
% Cleaned up and optimized 12/2019 by Aaron Levy to use only built-in
% matlab function (mostly replacing alphavol and inhull/inalpha with
% alphaShape and inShape and volume(alphaShape)
% Updated again late 2021,early 2022 ADL for use as PAINTshift of a whole field rather than
% xcorr of synapses

distance = distance.*psize;


edge=2; 
mv = 4;     
%   1/mv of the mean localization density will be used as cutoff 
%   to remove all high density peaks when calculation the best
%   overlapping.

%~~~Generate 3D distribution matrices for both channels

% Extract xyz coordinates of each matrix
x1 = A(:,1); 
y1 = A(:,2); 

x2 = B(:,1); 
y2 = B(:,2); 


% Shift xyz to have xyz=0 at min xyz position, then shift over by edge+1
x0 = min(min([x1;x2])); 
y0 = min(min([y1;y2])); 

x1 = x1 - x0 + edge + 1; 
y1 = y1 - y0 + edge + 1; 

x2 = x2 - x0 + edge + 1; 
y2 = y2 - y0 + edge + 1; 


% Convert xyz to renderpixel space
x1 = floor(x1 * psize / renderpixel);
y1 = floor(y1 * psize / renderpixel); 

x2 = floor(x2 * psize / renderpixel);
y2 = floor(y2 * psize / renderpixel); 


% Make arrays I1 and I2, which are pixelated versions of the synapse where
% the value at an index is the number of locs at that index
x0 = max(max([x1;x2])) + edge * psize / renderpixel; 
y0 = max(max([y1;y2])) + edge * psize / renderpixel; 

I1 = zeros(x0,y0);
for i = 1:length(x1)
    I1(x1(i),y1(i)) = I1(x1(i),y1(i)) + 1;
end
I2 = zeros(x0,y0);
for i = 1:length(x2)
    I2(x2(i),y2(i)) = I1(x2(i),y2(i)) + 1;
end

%~~~Generate mask for each channel
% Generate a mask the size of I1, set all points less than the synaptic
% area to nan to effectively crop the search radius of inShape, reset the
% mask values to 0, find the points of the mask that are inShape, and then
% populate the mask when the point is in shape. Finally reset all the nan
% to 0 for the full size mask.
mask1 = I1;
lim = min([x1;y1])-1;                 % this is the min value for xyz set above with some safety
mask1(1:lim,:) = nan;                  % set points that are below the minimum xyz positions to nan
mask1(:,1:lim) = nan;

mask1(~isnan(mask1)) = 1;
mask1(isnan(mask1)) = 0;

mask2 = I2;
lim = min([x2;y2])-1;                 % this is the min value for xyz set above with some safety
mask2(1:lim,:) = nan;                  % set points that are below the minimum xyz positions to nan
mask2(:,1:lim) = nan;

mask2(~isnan(mask2)) = 1;
mask2(isnan(mask2)) = 0;

%~~~Do the actual cross correlation

N1 = sum(sum(sum(I1.*mask1)));  % number of particles within channel 1
N2 = sum(sum(sum(I2.*mask2)));  % number of particles within channel 2
A1 = sum(sum(sum(mask1)));      % area of mask1
A2 = sum(sum(sum(mask2)));      % area of mask2
I1 = double(I1);   
I2 = double(I2);                % convert to double

NP = real(fftshift(ifftn(fftn(mask1,size(I1)+rmax).*conj(fftn(mask2,size(I1)+rmax))))); % Normalization for correct boundary conditions
tt = find(NP == 0); 
if size(tt,1) > 0 
    NP(tt) = -1e-15;            % if any values are zero change them to a very small value<0
end

G1 = A1*A2/N1/N2*real(fftshift(ifftn(fftn(I1.*mask1,size(I1)+rmax).*conj(fftn(I2.*mask2,size(I1)+rmax)))))./NP; % 2D G(r) with proper normalization
L0 = size(G1);
L = floor(L0 / 2 + 1);
NPmask = max(0,max(max(max(NP))) / 4-NP);
G1 = G1 + 0.0001;
G1(NPmask ~= 0) = 0;            % set non-relevant points to 0

% find xyzshift
if isempty(xyzshift)
    
    % remove high density peaks from the pixelated 3d loc distribution
    n = 11; 
    k = ones(n,n)./n^2; 
    I1a = convn(I1,k,'same');       % do n-dimensional convolution of k with I1 and I2
    I2a = convn(I2,k,'same'); 
    vA = (max(x1)-min(x1)) * (max(y1)-min(y1));
    d1 = length(x1) / vA;           % get average density of locs
    vB = (max(x2)-min(x2)) * (max(y2)-min(y2));
    d2 = length(x2) / vB;
    I1b = min(d1 / mv / (psize^2) * (renderpixel^2), I1a);  % remove high density peaks from I1 and I2
    I2b = min(d2 / mv / (psize^2) * (renderpixel^2), I2a);
    
    % convolve the mask, then calculate cross correlation
    n = 5; 
    k = ones(n,n)./n^2; 
    mask3 = convn(mask1,k,'same');
    mask4 = convn(mask2,k,'same'); 
    G0 = real(fftshift(ifftn(fftn(I1b.*mask3,size(I1b)+rmax).*conj(fftn(I2b.*mask4,size(I1b)+rmax))))); 
    
    % detect peak of correlation only in area with [distance] of center
    md = distance(2) / renderpixel; 
    G = G0 * 0; 
    G(L(1) - md:L(1) + md, L(2) - md:L(2) + md) = 1;
    [x, y] = ind2sub(size(G), find(G));
    d = pdist2([x y],L); 
    pG = find(G);
    G(pG(d > md)) = 0;
    G(pG(d < distance(1) / renderpixel)) = 0;

    G = G .* G0; 
    
    % get the xyz shift for max overlap
    Gm = max(max(max(G)));
    [mx, my] = ind2sub(size(G), find(G == Gm));         %find the peak of correlation
    xyzshift = [mx(1) - L(1), my(1) - L(2)]; 
    xyzshift = xyzshift .* (renderpixel / psize);
    
end
