%% function col = getColumns(headerstring)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% This function will return a struct col with the column indices for
% localization files loaded with Omniloader. For example, 
%   [headerstring,loc] = Omniloader('all_localizations_found0.txt');
%   col = getColumns(headerstring)
%col = 
%      frame: 8
%          x: 1
%          y: 2
%          z: [1×0 double]
%     sigmax: 3
%     sigmay: 6
%       phot: 10
%      error: 13
%        lpx: [1×0 double]
%        lpy: [1×0 double]
%        amp: 4
%         R2: 14
% This should work for matlab .txt files, ThunderSTORM MLE and Phasor 2d/3d
% .csv files, Picasso ThunderSTORM-formatted .csv files, and Picasso .hdf5
% files (2d/3d).
% 4/8/20 Aaron Levy
function col = getColumns(headerstring)

if isempty(headerstring)
    warning('Headerstring not loaded, please try again')
    return
end

headers = strsplit(headerstring,',');
col.frame = find(~cellfun('isempty',regexp(headers,'^frame','match')));
col.x = find(~cellfun('isempty',regexp(headers,'^x\s(\[)(px|nm)?(\])','match')));
col.y = find(~cellfun('isempty',regexp(headers,'^y\s(\[)(px|nm)?(\])','match')));
col.sigmax = find(~cellfun('isempty',regexp(headers,'^s(igma)?[x1]?\s(\[)(px|nm)(\])','match')));
col.sigmay = find(~cellfun('isempty',regexp(headers,'^s(igma)?[y2]\s(\[)(px|nm)(\])','match')));
col.error = find(~cellfun('isempty',regexp(headers,'^(error|uncertainty_xy)?\s(\[)(px|nm)(\])','match')));
col.phot = find(~cellfun('isempty',regexp(headers,'^(photons|intensity)?(\s[photon(s)?)])?','match')));
col.lpx = find(~cellfun('isempty',regexp(headers,'^lpx\s(\[)(px|nm)?(\])','match')));
col.lpy = find(~cellfun('isempty',regexp(headers,'^lpy\s(\[)(px|nm)?(\])','match')));
col.amp = find(~cellfun('isempty',regexp(headers,'^amplitude','match')));
col.R2 = find(~cellfun('isempty',regexp(headers,'^R\^2','match')));
col.channel = find(~cellfun('isempty',regexpi(headers,'^channel','match')));
col.bg = find(~cellfun('isempty',regexpi(headers,'^bg','match')));
col.ellip = find(~cellfun('isempty',regexpi(headers,'^ellipticity','match')));
col.netgrad = find(~cellfun('isempty',regexpi(headers,'^net_gradient','match')));
col.len = find(~cellfun('isempty',regexpi(headers,'^len','match')));
col.n = find(~cellfun('isempty',regexpi(headers,'^n$','match')));
col.photrate = find(~cellfun('isempty',regexpi(headers,'^photon_rate','match')));
col.dz = find(~cellfun('isempty',regexpi(headers,'^d_zcalib','match')));

col.groups = find(~cellfun('isempty',regexpi(headers,'^group(s)?','match')));
col.convex_hull = find(~cellfun('isempty',regexpi(headers,'^convex_hull','match')));
col.area = find(~cellfun('isempty',regexpi(headers,'^area','match')));
col.com_x = find(~cellfun('isempty',regexpi(headers,'^com_x','match')));
col.com_y = find(~cellfun('isempty',regexpi(headers,'^com_y','match')));
col.std_frame = find(~cellfun('isempty',regexpi(headers,'^std_frame','match')));
col.mean_frame = find(~cellfun('isempty',regexpi(headers,'^mean_frame','match')));
col.std_x = find(~cellfun('isempty',regexpi(headers,'^std_x','match')));
col.std_y = find(~cellfun('isempty',regexpi(headers,'^std_y','match')));
col.LLrel = find(~cellfun('isempty',regexpi(headers,'^LLrel','match')));

col.synnum = find(~cellfun('isempty',regexpi(headers,'^synnum','match')));


% If you loaded a Phasor file, then sigma is not useful since it's not
% sigma from a fit. The phasor files have 'sigma1' labeles, but do not have
% an error column, so in that case set the x and y sigma to empty.
if ~isempty(regexp(headers{col.sigmax},'^(sigma1)*','match')) && isempty(col.error) && isempty(col.lpx)
    col.sigmax = double.empty(1,0); 
    col.sigmay = double.empty(1,0);  
end

fn = fieldnames(col);
tf = cellfun(@(c) isempty(col.(c)), fn);
col = rmfield(col, fn(tf));

col.z = find(~cellfun('isempty',regexp(headers,'^z\s(\[)(px|nm)?(\])','match')));


end
