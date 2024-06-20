%% function vals = extractNCAnalyses(s, dataset, datatype, field, varargin)
% function to pull out variables from NanostructureAnalysisSuitePD.
% INPUTS
% s = struct output of the code (aka AllData)
% dataset = dataset within s to analyze (aka 'Kaeser')
% datatype = 'Syn' or 'NC'
% field = fieldname to extract like 'Channel_3_NCnum'
% OUTPUTS
% vals = array of data from field
% 5/1/23 ADL
function vals = extractNCAnalyses(s, dataset, datatype, field, varargin)

p = inputParser;
addRequired(p,'s',@isstruct)
addRequired(p,'dataset',@ischar)
opts = {'Syn','NC'};
addRequired(p,'datatype',@(x)mustBeMember(x,opts))
addRequired(p,'field',@ischar)
addOptional(p,'POI','',@ischar)
parse(p,s,dataset,datatype,field,varargin{:})
POI = p.Results.POI;

% s = AllData;
% dataset = 'Kaeser';
% datatype = 'Syn'; % or 'NC'
% field = 'Channel_3_NCnum';

subs = s.([dataset '_' datatype '_Dataset']);

if strcmp(datatype,'NC')
    subs = subs.([POI '_Data']);
end

if sum(ismember(fieldnames(subs),field)) == 1
    vals = cell2mat({subs(:).(field)}');
else
    error('Unavailable fieldname, try again')
end

end
