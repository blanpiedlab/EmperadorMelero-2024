%% dilatedCluster = dilateCluster(clusterLoc,ringwidth,psize,varargin)
% This function takes a synaptic cluster and dilates it by distance
% ringwidth to expand the cluster (or contract if ringwidth<0)
% REQUIRED INPUTS
% clusterLoc = m x 2 array of xy positions of cluster to dilate
% ringwidth = distance you want to dilate in nm; negative values contract
% psize = image pixelsize, ie 160
% OPTIONAL NAME/VALUE INPUTS
% figflag = default 0; set to 1 to display plot of OG and dilated cluster
% OUTPUTS
% dilatedCluster = polyshape of new, dilated cluster
% Written 4/2023 ADL; dilation code from PD
function dilatedCluster = dilateCluster(clusterLoc,ringwidth,psize,varargin)
    
    p = inputParser;
    addRequired(p,'clusterLoc',@isnumeric)
    addRequired(p,'ringwidth',@isnumeric)
    addRequired(p,'psize',@isnumeric)
    addParameter(p,'figflag',0,@isnumeric)
    parse(p,clusterLoc,ringwidth,psize,varargin{:})
    
    if size(clusterLoc,2)>2
        error('Load just xy positions for clusterLoc rather than entire loc table');
    end
    
    % create a dilated version of the cluster
    testAlphaShape = alphaShape(p.Results.clusterLoc(:,1), p.Results.clusterLoc(:,2),2);
    [~,synboundary] = boundaryFacets(testAlphaShape);
    testPolyShape = polyshape(synboundary(:,1), synboundary(:,2));
    ringradiusnmtopx = p.Results.ringwidth/p.Results.psize;
    dilatedCluster = polybuffer(testPolyShape, ringradiusnmtopx, 'JointType','round');
  
    if p.Results.figflag == 1
        
        figure
        hold on
        plot(testPolyShape,'EdgeAlpha',0)
        plot(dilatedCluster,'FaceAlpha',0.2)
        legend({'Original Cluster','Dilated Cluster'})
        
    end

end