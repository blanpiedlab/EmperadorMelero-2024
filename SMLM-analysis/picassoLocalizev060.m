%% function [outpath,yamlpath] = picassoLocalize(impath,boxsize,fitmethod,gradient,drift,camera,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% This function will run picasso localize through the windows cmd line
% Must have picassopip env setup on the windows path for this to work.
% If loading ome.tif, picasso will grab all the same name ome.tif from the
% folder and localize together (ie output off PAINT scope). If coming from
% the STORM scope, save out as raw first and feed this the raw file. YAML
% must also exist.
% REQUIRED INPUTS:
% - impath (char): full path to where the tif or raw file is
%    ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A.raw'   
% - psize (num): camera pixel size in nm/px (ie 160)
% - boxsize (num): pixel size of box for fitting (ie 7 for 2d and 11 for 3d)
% - fitmethod (char): choose from {'mle','lq','lq-gpu','lq-3d','lq-gpu-3d','avg'}
% - gradient (num): min net gradient for filtering spots
% - drift (1 or 0): set to 0 to turn off automatic drift correction after localization
% - camera (char): enter 'ixon','prime-16-fullwell','prime-12-fullwell',...
%         'prime-12-balanced',or 'prime-12-sensitivity' to have those camera params loaded
%     optionally enter 'custom' for camera and then enter 4 arguments, in order,
%         for baseline, sensitivity, EM gain, quantumn efficiency, and pixel size in
%         that order (as numbers)
%             ie ...'custom',100,0.5,300,0.95,160
%     example where impath = 'Z:\Aaron\testims\testim.tif'
%       picassoLocalize(impath,7,'lq-gpu',5000,0,'prime-12-balanced');
%     example: picassoLocalize(impath,7,'lq-gpu',5000,0,'custom',100,0.5,300,0.75);
% NAME/VALUE OPTIONAL INPUTS:
% - verbose: default 1, set to 0 to suppress cmd line output from Python
% - magfac: set to value between 0>x>1 to define z fitting magnification
% factor for 3d images
% - zcalipath: full path to z calibration yaml file
% - suffix: add a suffix to the file for example '_right' will yield
%       'filename_right_locs.hdf5'
% - roi: list an roi in a string format minx,miny,maxx,maxy (ie top left
%   and bottom right coords) so '0,0,256,300' would do the left side and
%   '256,0,512,300' would do the right side.
% OUTPUTS
% - outpath: full path to '_locs.hdf5' file
% - yamlpath: full path to '_locs.yaml' file

function [outpath,yamlpath] = picassoLocalizev060(impath,boxsize,fitmethod,gradient,drift,camera,varargin)

    % Check for picasso command to appear in cmd line
    [status,~] = system('picasso localize -h');
    if status ~= 0
        error('Something is not right in the system environment variable for cmd to find the path to picassopip env.')
    end

    % Parse some inputs
    p = inputParser;
  
    addRequired(p,'impath',@ischar)
    addRequired(p,'boxsize',@isnumeric)
    opts = {'mle','lq','lq-gpu','lq-3d','lq-gpu-3d','avg'};
    addRequired(p,'fitmethod',@(x)mustBeMember(x,opts))
    addRequired(p,'gradient',@isnumeric)
    addRequired(p,'drift',@isnumeric)
    addRequired(p,'camera',@ischar)
    addOptional(p,'bl',@isnumeric)
    addOptional(p,'s',@isnumeric)
    addOptional(p,'ga',@isnumeric)
    addOptional(p,'qe',@isnumeric)
    addOptional(p,'px',@isnumeric)
    addParameter(p,'verbose',1,@isnumeric)
    fun = @(x) (x >= 0) && (x <= 1);
    addParameter(p,'magfac',nan,fun)
    addParameter(p,'zcalipath','',@isfile)
    addParameter(p,'roi','',@ischar)
    addParameter(p,'suffix','',@ischar)
    parse(p,impath,boxsize,fitmethod,gradient,drift,camera,varargin{:})
    
    % Check input file is tif or raw
    [folder,file,ext] = fileparts(impath);
    if sum(strcmp(ext,{'.tif','.raw','.nd2'}))==0
        error('You did not feed this function a tif, raw, or nd2 file, try again.')
    end

    impath = p.Results.impath;
    b = num2str(p.Results.boxsize);
    g = num2str(p.Results.gradient);
    d = num2str(p.Results.drift);
    a = p.Results.fitmethod;
    mf = num2str(p.Results.magfac);
    zc = p.Results.zcalipath;
    roi = p.Results.roi;
    sfx = p.Results.suffix;
    
    % Get the camera settings
    switch p.Results.camera 
        case 'ixon'
            bl = '100';
            s = '0.12';
            ga = '1';
            qe = '0.95';
            px = '160';
        case 'prime-16-fullwell'
            bl = '100';
            s = '0.8';
            ga = '1';
            qe = '0.95';
            px = '100';
        case 'prime-12-fullwell'
            bl = '100';
            s = '2.56';
            ga = '1';
            qe = '0.95';
            px = '100';
        case 'prime-12-balanced'
            bl = '100';
            s = '1.16';
            ga = '1';
            qe = '0.95';
            px = '100';
        case 'prime-12-sensitivity'
            bl = '100';
            s = '0.57';
            ga = '1';
            qe = '0.95';
            px = '100';
        case 'custom'
            bl = num2str(varargin{1});
            s = num2str(varargin{2});
            ga = num2str(varargin{3});
            qe = num2str(varargin{4});
            px = num2str(varargin{5});
        otherwise
            error('You probably input the camera name wrong, try again')
    end
    
    % Do the actual localization
    if ~contains(fitmethod,'3d')
        
        S = ['picasso localize "' impath '" -px ' px ' -b ' b ' -a ' a ...
                ' -g ' g ' -d ' d ' -bl ' bl ' -s ' s ' -ga ' ga ' -qe ' qe];
        
    elseif contains(fitmethod,'3d')
        
        S = ['picasso localize "' impath '" -px ' px ' -b ' b ' -a ' a ...
            ' -g ' g ' -d ' d ' -bl ' bl ' -s ' s ' -ga ' ga ' -qe ' qe...
            ' -mf ' mf ' -zc "' zc '"'];
        
    end
    
    if ~isempty(roi)
        
        S = [S ' -r ' roi];
        
    end
    
    if ~isempty(sfx)
        
        S = [S ' -sf ' sfx];
        
    end

    out = 'Localizing with pixelsize %s box size %s, min gradient %s, fit method %s, and camera %s.\n';
    fprintf(out,px,b,g,a,p.Results.camera)

    switch p.Results.verbose
        case 1
            [status,~] = system(S,'-echo'); % do localization and print python cmd
        case 0
            [status,~] = system(S); % do localization without print
    end
    
    % Output whether it worked
    if status == 0
        fprintf('Localization complete!\n')
    else
        fprintf('There was an error in localization, try again!\n')
    end

    % output file
    if ~isempty(sfx)
        outpath = fullfile(folder,[file sfx '_locs.hdf5']);
        yamlpath = fullfile(folder,[file sfx '_locs.yaml']);
    else
        outpath = fullfile(folder,[file '_locs.hdf5']);
        yamlpath = fullfile(folder,[file '_locs.yaml']);
    end

end
