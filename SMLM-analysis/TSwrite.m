%% TSwrite(loc,savename,header,varargin)
% TSwrite is a catchall function for saving out ThunderSTORM-style
% localization files from Matlab. The function will take the localization
% table you hae made in Matlab as well as a savename, and a header (either
% one you provide or one pulled from a previous file), and then use fprintf
% to save out the new file with a default precision of %.5f (same as
% default ThunderSTORM output).
%
% REQUIRED INPUTS
% *loc: The loc table that you want to export, ie nloc x 10 for phasor, etc
% *savename: The filename you want this to save under. Must include csv, ie
%   'TS_saved.csv'
% *header: The header line source. Can take one of two forms:
%       1: enter a previous CSV filename from the same directory
%       (ie, 'TS_output.csv') and TSwrite will automatically pull its
%       header and use it here.
%       2: enter a string of header titles, with columns separated by
%       commas without a space, as below for example for Phasor
%       header = 'id,frame,x [nm],y [nm],z [nm],intensity [photon],offset [photon],bkgstd [photon],sigma1 [nm],sigma2 [nm]'
% OPTIONAL INPUT
% *name/value pairs:
%     -'fprec', precision: default precision is 5 (ie, %.5f), but can be
%     changed with this. 
%     -'convertToPixels',pixelsize: Will divide any columns with [nm] by
%     pixelsize to convert to camera pixel space and save in pixel space
%     -'convertToNm',pixelsize: will multiply any columns with [px] by
%     pixelsize to convert from camera pixel space to nm and save in nm.
%     -'verbose', 0: default is 1 (on). set to 0 to suppress all outputs
%     besides the saved column order.
% OUTPUT
% loc will be saved to cd and called savename, with precision fprec
% (default 5) and with header defined by header variable

% Version history
% v1.0 Aaron Levy 8/30/19
% v2.0 Aaron Levy 3/26/20, added name/value pairs and updated some text and
% inputs

% % test variables
% autoheaderfile = 'TS_float5.csv';
% xcol=3;
% ycol=4;
% photcol=6;
% filter=0;
% loc = TSread(autoheaderfile,xcol,ycol,photcol,filter);
% savename = 'TS_float5_saved.csv';
% headerstring = 'id,frame,x [nm],y [nm],z [nm],intensity [photon],offset [photon],bkgstd [photon],sigma1 [nm],sigma2 [nm]';
% fprec = 5;

function TSwrite(loc,savename,header,varargin)

p = inputParser;
addRequired(p,'loc');
addRequired(p,'savename',@ischar);
addRequired(p,'header',@ischar);
addParameter(p,'fprec',5);
addParameter(p,'convertToPixels',0);
addParameter(p,'convertToNm',0);
addParameter(p,'verbose',1);
parse(p,loc,savename,header,varargin{:});

loc = p.Results.loc;
savename = p.Results.savename;
header = p.Results.header;
fprec = p.Results.fprec;
convertToPixels = p.Results.convertToPixels;
convertToNm = p.Results.convertToNm;
verbose = p.Results.verbose;

if isempty(header)
    
    warning('You need to specify either an autoheader filename or the actual header as a string');
    return;
    
elseif endsWith(header,{'.csv','.txt','.hdf5'})
   
    % auto get header from the original file
    S = 'Saving locs using header from provided filename...\n';
    if verbose; fprintf(S); end
    
    filename = header;
    
    [header, ~] = Omniloader(filename);
   
else

    S = 'Saving locs using header provided in function input...\n';
    if verbose; fprintf(S); end

end

% exit the function if the header doesn't have the same number of
% columns as the loc variable to save does
test = length(strsplit(header,','));
ncol = size(loc,2);
if test ~= ncol
    
    warning(['Your header does not have the right number of columns,' ...
    ' correct the header and try again']);
    return;
    
end

% Optionally convert columns labeled [nm] to camera pixels
if convertToPixels ~= 0
    
    headers = strsplit(header,',');
    
    xyzcols = find(~cellfun('isempty',regexp(headers,'^(x|y|z)?\s(\[)nm(\])','match')));
    sigmacols = find(~cellfun('isempty',regexp(headers,'^s(igma)?(x|y|\s)?\s(\[)nm(\])','match')));
    cols = [xyzcols sigmacols];
    for i = 1:length(cols)
        thiscol = cols(i);
        loc(:,thiscol) = loc(:,thiscol)./convertToPixels;
        headers(thiscol) = replace(headers(thiscol),'[nm]','[px]');
    end
    
    header = strjoin(headers,',');

    S = ['\nYou have converted xyz and sigma that were in nm to camera pixels '...
        'with pixel size %d nm.\n ***NOTE*** that error is still in nm, check headerstring.\n'];
    if verbose; fprintf(S,convertToPixels); end

    
end

% Optionally convert columns labeled [px] to nm

if convertToNm ~= 0
    
    headers = strsplit(header,',');
    
    xyzcols = find(~cellfun('isempty',regexp(headers,'^(x|y|z)?\s(\[)px(\])','match')));
    sigmacols = find(~cellfun('isempty',regexp(headers,'^s(igma)?(x|y|\s)?\s(\[)px(\])','match')));
    errorcols = find(~cellfun('isempty',regexp(headers,'^lp(x|y)?\s(\[)px(\])','match')));
    cols = [xyzcols sigmacols errorcols];
    for i = 1:length(cols)
        thiscol = cols(i);
        loc(:,thiscol) = loc(:,thiscol).*convertToNm;
        headers(thiscol) = replace(headers(thiscol),'[px]','[nm]');
    end
    
    header = strjoin(headers,',');
    
    S = ['\nYou have converted xyz and sigma that were in camera pixels '...
        'to nm, with pixel size %d nm.\n***NOTE***Error should also be in nm.\n'];
    if verbose; fprintf(S,convertToNm); end
    
end

% Display the saved column names
cnames = strsplit(header,',');
nums = strsplit(num2str(1:length(cnames)),' ');
shortstring1 = '\t';
shortstring2 = ': %s\n';
string = strjoin(strcat(shortstring1,nums,shortstring2));
fprintf('Your saved columns are:\n');
fprintf(string,cnames{:});

% generate the string of modifiers based on number of columns
shortstring = ['%.' num2str(fprec) 'f'];
string = strjoin(repmat({shortstring},1,ncol));
string = insertBefore(string,' ',',');

% save out the data
tic
fid = fopen(savename,'wt');
fprintf(fid,[header '\n']);
fprintf(fid,[string '\n'],loc.');
fclose(fid);

S = 'Saved %.0f locs in %.1f seconds with precision %.0f.\n';
if verbose; fprintf(S,length(loc),toc,fprec); end

end %function
