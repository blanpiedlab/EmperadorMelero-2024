function randcluster = get_cluster_randomized_ADL_PD2D_specifyShape(xyz,factor,shp,flag)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy, Poorna Dharmasri, Aihui Tang
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% function used to randomize all localizations inside a cluster A
% input:
% A: [x,y,z]
% factor: density adjusting factor.
% pisze = camera pixel size (nm/px)
% flag 1 for figure

% output:
% randcluster: new [x,y,z] with all localizations randomized



if nargin ~= 4
    flag = 0;
end

% Extract xyz, number of locs, and alphaShape volume
x = xyz(:,1); 
y = xyz(:,2); 

n = size(x,1); 
%shp = alphaShape(xyz,150/psize);
v = area(shp);

minshapex = min(shp.Points(:,1));
maxshapex = max(shp.Points(:,1));
minshapey = min(shp.Points(:,2));
maxshapey = max(shp.Points(:,2));

% Generate random localizations at the same loc density as the real synapse
% inside a rectangle bounding the synapse 
md = n*factor/v;                                            % mean loc density in synapse
mv = (maxshapex-minshapex) * (maxshapey-minshapey);   % mean rectangular volume
num = md*mv;                                                % number of locs needed for synaptic loc density in the rectangle
num = floor(sqrt(num)*3);                                   % honeslty, not sure what this factor is. some kind of dimensional reduction
r0 = random('unif',minshapex,maxshapex,num);     randx = r0(:);   % generate a num x num array of random numbers between min(x) and max(x)
r0 = random('unif',minshapey,maxshapey,num);     randy = r0(:);

randxyz=[randx'; randy']';                          % create a (num*num) x 3 array of random values spanning the rectangle fitting the synapse volume

% Crop the region to the synapse volume and the correct number of points
tf = inShape(shp,randxyz(:,1),randxyz(:,2));   % find which random points are in the alphaShape
randxyz = randxyz(tf == 1,:);                               % reset randxyz to just points in the alphaShape
randcluster = randxyz(1:floor(n*factor),:);                 % takes the first n*factor (same # as synaptic locs modified by factor) random values

if flag == 1
    tiledlayout(1,2)
    nexttile
    scatter(x,y,'.b')
    title('original cluster')
    axis equal
    nexttile
    scatter(randcluster(:,1),randcluster(:,2),'.b');
    title('random cluster');
    axis equal
end
