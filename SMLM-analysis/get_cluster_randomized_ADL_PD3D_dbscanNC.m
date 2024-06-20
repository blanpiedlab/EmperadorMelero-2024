function randcluster = get_cluster_randomized_ADL_PD3D_dbscanNC(xyz,factor,psize,flag)
% function used to randomize all localizations inside a cluster A
% input:
% A: [x,y,z]
% factor: density adjusting factor.
% pisze = camera pixel size (nm/px)
% flag 1 for figure

% output:
% randcluster: new [x,y,z] with all localizations randomized

% Updated 12/2019 by ADL to use all builtin matlab functions (ie alphaShape
% instead of alphaVol) and cleand up a bit.

if nargin ~= 3
    flag = 0;
end

% Extract xyz, number of locs, and alphaShape volume
x = xyz(:,1); 
y = xyz(:,2); 
z = xyz(:,3);

n = size(x,1); 
shp = alphaShape(x,y,z,150/psize);
v = volume(shp);

% Generate random localizations at the same loc density as the real synapse
% inside a rectangle bounding the synapse 
md = n*factor/v;                                            % mean loc density in synapse
mv = (max(x)-min(x)) * (max(y)-min(y)) * (max(z)-min(z));   % mean rectangular volume
num = md*mv;                                                % number of locs needed for synaptic loc density in the rectangle
num = floor(sqrt(num)*3);                                   % honeslty, not sure what this factor is. some kind of dimensional reduction
r0 = random('unif',min(x),max(x),num);     randx = r0(:);   % generate a num x num array of random numbers between min(x) and max(x)
r0 = random('unif',min(y),max(y),num);     randy = r0(:);
r0 = random('unif',min(z),max(z),num);     randz = r0(:);
randxyz=[randx'; randy';randz']';                          % create a (num*num) x 3 array of random values spanning the rectangle fitting the synapse volume

% Crop the region to the synapse volume and the correct number of points
tf = inShape(shp,randxyz(:,1),randxyz(:,2),randxyz(:,3));   % find which random points are in the alphaShape
randxyz = randxyz(tf == 1,:);                               % reset randxyz to just points in the alphaShape
randcluster = randxyz(1:floor(n*factor),:);                 % takes the first n*factor (same # as synaptic locs modified by factor) random values

if flag == 1
    scatter3(randcluster(:,1),randcluster(:,2));
    title('random cluster');
end