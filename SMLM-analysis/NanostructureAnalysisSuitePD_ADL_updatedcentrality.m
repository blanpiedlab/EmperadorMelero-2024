% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Poorna Dharmasri
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Note that this code outputs a lot of analyses, not all of which were used in the paper above.
%%
clear
clc
close all
tic
warning('off','MATLAB:alphaShape:DupPointsBasicWarnId');

%%%%%%%%%%% ENTER INFO IN THIS SECTION BEFORE RUNNING SCRIPT %%%%%%%%%%%%%

%%%% DONT FORGET TO SET THIS UP AHEAD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ensure that NCparams are saved by protein name in the dataset folder in an nx3 array
%%% where NCparams(:,1) = synnum; NCparams(:,2) = epsilon; NCparams(:3) = minpts
%%% (ie 'PSD95_NCparams.mat')
%%%% DONT FORGET TO SET THIS UP AHEAD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Name the channels in the ProteinsInQuestionIDs cell array.
%%% Set the channels you want to analyze NCs in by modifying the
%%% iq_channels array with the channel numbers separated by a space. Set
%%% the channels whose distributions you want to measure against an NC in
%%% question using the target protein or tp_channels array. For example,
%%% where an iq_channel == tp_channel, you get auto enrichment. Where it
%%% doesn't equal, you get cross enrichments.

ProteinsInQuestionIDs = {'Bsn', 'PSD95', 'Munc', 'Cav'};
TargetProteinIDs = ProteinsInQuestionIDs;
iq_channels = [1 2 3 4];
tp_channels = [1 2 3 4];
relevant_channels = unique([iq_channels, tp_channels]);

shapeinstructionarray = [1 1; 2 2; 3 1; 4 1]; %This array  is used to instruct which
% protein's shape to use to randomize within for analytical reasons. The
% first column indicates the protein of interest or target protein, depending on the analysis,
% and the second column instructs which protein's shape to use.

%%% Give some info on which columns hold the folloiwng information, where
%%% channelcol is the column designating which protein species is imaged,
%%% synapsecol is the synapse number column
%%% xcol and ycol are x and y coordinate respectively
channelcol = 19;
synapsecol = 20;
xcol = 5;
ycol = 6;
psize = 160;
filename = 'synapse*.txt'; %Sets the base file name for each loc file. Change as necessary
randflag = 0; %IF 1, do analyses with randomized NCs to compare to real data. IF 0, don't do rand NCs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('WELCOME TO THE NANOSTRUCTURE ANALYSIS SUITE WE HOPE YOU ENJOY YOUR STAY\n\n\n')
%Create the megastructure that will house all the data, and start setting
%up folder searching for files we want to analyze
AllData = struct;

% Set up some variables for the analyses
randruns = 100;
radius = (5:10:330)/psize; % center of bins to use for NC enrichment see below for conversion to edges
d = diff(radius)/2;
radius = [radius(1)-d(1), radius(1:end-1)+d, radius(end)+d(end)];
radius(2:end) = radius(2:end)+ eps(radius(2:end));

% Identify the dataset folders to loop through
maindir = pwd;
ds = dir;
dfolders = {ds([ds.isdir]).name};
dfolders = dfolders(~ismember(dfolders,{'.','..'}));
datasetdir = cell(numel(dfolders),1);
for fn = 1:size(dfolders,1)
    datasetdir{fn} = fullfile(maindir,dfolders{fn});
end

%Start the meta loop that loops through all the folders
for foldnum = 1:length(datasetdir)
    %Because this was designed for a project where we had multiple
    %datasets,(e.g. dataset 1, dataset 2, etc), we run through
    %each of those and create a structure to house the data from each of
    %those folders.
    thispath = datasetdir{foldnum};
    [~,ds_ID,~] = fileparts(thispath);
    fprintf('Now loading data from the %s dataset...\n',ds_ID)
    
    ThisDataSetNCData = struct;
    ThisDataSetSynData = struct;
    searchterm = fullfile(thispath,filename);
    files = dir(searchterm);
    
    %Following section concatenates all the synapse.txt files in a given
    %folder in order, and renumbers them from 1 to n, where n is all the synapses in
    %that dataset (e.g. 149 for Kaeser data)
    allsyns = [];
    for filenum =  1:size(files,1)
        
        thisfilename = ['synapse' num2str(filenum) '.txt'];
        tempfile = load(fullfile(thispath,thisfilename));  %use for .txt files
        
        if filenum == 1
            allsyns = [allsyns; tempfile];
        elseif filenum >1
            maxsyn = max(allsyns(:,synapsecol));
            tempfile(:,synapsecol) = tempfile(:,synapsecol) + maxsyn;
            allsyns = [allsyns; tempfile];
        end
        
    end
    
    % set up some additional column variables for peak and ncID storage
    % based on the size of allsyns
    [nrow,ncol] = size(allsyns);
    ncIDcol = ncol+1;
    peakxcol = ncol+2;
    peakycol = ncol+3;
    
    % Load the NC params that were precalculated
    % NOTE OUTLIERS assumes that your area cutoffs were given in nm^2
    % rather than pixels^2. otherwise adjust.
    NCparams = struct;
    tempoutliers = struct2cell(load(fullfile(thispath,'NCareathresh.mat')));
    tempoutliers = tempoutliers{1};
    outliers = struct;
    for i = 1:numel(relevant_channels)
        
        thisch = ProteinsInQuestionIDs{relevant_channels(i)};
        filename = fullfile(thispath,[thisch '_NCparams.mat']);
        temp = load(filename);
        temp = struct2cell(temp);
        NCparams.(['Channel_' num2str(relevant_channels(i))]) = temp{1};
        outliers.(['Channel_' num2str(relevant_channels(i))]) = tempoutliers.(thisch)./psize^2;
    end
    
    %Below, create structure for each channel (i.e. Protein) based on if we
    %want to analyze it or not for this dataset.
    
    if sum(ismember(iq_channels,1))
        Channel_1_struct = struct;
    end
    if sum(ismember(iq_channels,2))
        Channel_2_struct = struct;
    end
    if sum(ismember(iq_channels,3))
        Channel_3_struct = struct;
    end
    if sum(ismember(iq_channels,4))
        Channel_4_struct = struct;
    end
    %Initialize a count for the number of NCs analyzed for each channel
    %because each row of each structure will be an individual NC
    
    channel1nctotal = 1;
    channel2nctotal = 1;
    channel3nctotal = 1;
    channel4nctotal = 1;
    
    fprintf('\tLoading complete, begin analysis\n')
    outliersyns = struct2cell(load(fullfile(thispath,'outliersynapses.mat')));
    outliersyns = outliersyns{1};
   
    % Loop across synapses
    errorhandler = struct; % initialize for at least minimal error handling
    for synnum = 1:max(allsyns(:,synapsecol))
        
        if ismember(synnum, outliersyns)
            fprintf('We got ourselves an outlier here folks, skipping to the next synapse!\n')
            continue
        end
        
        try % we are try/catching the entire synapse loop just for some basic error handling to ID any issues
            
            fprintf('We are on synapse %d\n',synnum)
            fprintf('Analyzing synapse data...')
            %%%% First calculate all the synaptic information
            
            % Grab the locs for each channel for this synapse
            thissyn = allsyns(allsyns(:,synapsecol)==synnum,:);
            thisch1 = thissyn(thissyn(:,channelcol)==1,:);
            thisch2 = thissyn(thissyn(:,channelcol)==2,:);
            thisch3 = thissyn(thissyn(:,channelcol)==3,:);
            thisch4 = thissyn(thissyn(:,channelcol)==4,:);
            
            % Get autocorrelations per channel
            ch1AC = get_autocorr_3d_ADL_PD2D(thisch1(:,[xcol ycol]), 5, 50, psize, 0);
            ch2AC = get_autocorr_3d_ADL_PD2D(thisch2(:,[xcol ycol]), 5, 50, psize, 0);
            ch3AC = get_autocorr_3d_ADL_PD2D(thisch3(:,[xcol ycol]), 5, 50, psize, 0);
            ch4AC = get_autocorr_3d_ADL_PD2D(thisch4(:,[xcol ycol]), 5, 50, psize, 0);
            
            % Get all combinations of crosscorrelations
            crosscorr12 = get_crosscorr_2d_ADL(thisch1(:,[xcol ycol]), thisch2(:,[xcol ycol]),psize,5,50);
            crosscorr13 = get_crosscorr_2d_ADL(thisch1(:,[xcol ycol]), thisch3(:,[xcol ycol]),psize,5,50);
            crosscorr14 = get_crosscorr_2d_ADL(thisch1(:,[xcol ycol]), thisch4(:,[xcol ycol]),psize,5,50);
            crosscorr23 = get_crosscorr_2d_ADL(thisch2(:,[xcol ycol]), thisch3(:,[xcol ycol]),psize,5,50);
            crosscorr24 = get_crosscorr_2d_ADL(thisch2(:,[xcol ycol]), thisch4(:,[xcol ycol]),psize,5,50);
            crosscorr34 = get_crosscorr_2d_ADL(thisch3(:,[xcol ycol]), thisch4(:,[xcol ycol]),psize,5,50);
            
            % Get synapse areas (for Kaeser data some of this is not going to
            % be meaningful but calculate it anyway since we need to set up the
            % synapse shape storage cell array anyway
            synapseshapestorage = cell(4,1);
            synapseshapestorage{1} = alphaShape(thisch1(:,xcol),thisch1(:,ycol),'HoleThreshold',100000);
            synapseshapestorage{2} = alphaShape(thisch2(:,xcol),thisch2(:,ycol),'HoleThreshold',100000);
            synapseshapestorage{3} = alphaShape(thisch3(:,xcol),thisch3(:,ycol),'HoleThreshold',100000);
            synapseshapestorage{4} = alphaShape(thisch4(:,xcol),thisch4(:,ycol),'HoleThreshold',100000);
            
            ch1synarea = area(synapseshapestorage{1});
            ch2synarea = area(synapseshapestorage{2});
            ch3synarea = area(synapseshapestorage{3});
            ch4synarea = area(synapseshapestorage{4});
            
            % Get centrality measures
            [ch1centrality, ch1centralitysqrt, ch1watanabecentrality] = Centrality(thisch1(:,[xcol ycol]),thisch1(:,[xcol ycol]),'unitcircle',1);
            [ch2centrality, ch2centralitysqrt, ch2watanabecentrality] = Centrality(thisch2(:,[xcol ycol]),thisch2(:,[xcol ycol]),'unitcircle',1);
            [ch3centrality, ch3centralitysqrt, ch3watanabecentrality] = Centrality(thisch1(:,[xcol ycol]),thisch3(:,[xcol ycol]),'unitcircle',1);
            [ch4centrality, ch4centralitysqrt, ch4watanabecentrality] = Centrality(thisch1(:,[xcol ycol]),thisch4(:,[xcol ycol]),'unitcircle',1);
            
            % Dump all into struct
            ThisDataSetSynData(synnum).('SynNum') = synnum;
            ThisDataSetSynData(synnum).('Week') = thissyn(1,1);
            ThisDataSetSynData(synnum).('Day') = thissyn(1,2);
            ThisDataSetSynData(synnum).('Region') = thissyn(1,3);
            ThisDataSetSynData(synnum).('Channel_1_locs') = thisch1;
            ThisDataSetSynData(synnum).('Channel_2_locs') = thisch2;
            ThisDataSetSynData(synnum).('Channel_3_locs') = thisch3;
            ThisDataSetSynData(synnum).('Channel_4_locs') = thisch4;
            
            ThisDataSetSynData(synnum).('Channel_1_SynArea') = ch1synarea;
            ThisDataSetSynData(synnum).('Channel_1_LocNum') = size(thisch1,1);
            ThisDataSetSynData(synnum).('Channel_1_SynDensity') = size(thisch1,1)/ch1synarea;
            ThisDataSetSynData(synnum).('Channel_1_AutoCorr') = ch1AC;
            ThisDataSetSynData(synnum).('Channel_1_Centrality') = ch1centrality;
            ThisDataSetSynData(synnum).('Channel_1_CentralitySqrt') = ch1centralitysqrt;
            ThisDataSetSynData(synnum).('Channel_1_WatanabeCentrality') = ch1watanabecentrality;
            
            ThisDataSetSynData(synnum).('Channel_2_SynArea') = ch2synarea;
            ThisDataSetSynData(synnum).('Channel_2_LocNum') = size(thisch2,1);
            ThisDataSetSynData(synnum).('Channel_2_SynDensity') = size(thisch2,1)/ch2synarea;
            ThisDataSetSynData(synnum).('Channel_2_AutoCorr') = ch2AC;
            ThisDataSetSynData(synnum).('Channel_2_Centrality') = ch2centrality;
            ThisDataSetSynData(synnum).('Channel_2_CentralitySqrt') = ch2centralitysqrt;
            ThisDataSetSynData(synnum).('Channel_2_WatanabeCentrality') = ch2watanabecentrality;
            
            ThisDataSetSynData(synnum).('Channel_3_SynArea') = ch3synarea;
            ThisDataSetSynData(synnum).('Channel_3_LocNum') = size(thisch3,1);
            ThisDataSetSynData(synnum).('Channel_3_SynDensity') = size(thisch3,1)/ch3synarea;
            ThisDataSetSynData(synnum).('Channel_3_AutoCorr') = ch3AC;
            ThisDataSetSynData(synnum).('Channel_3_Centrality') = ch3centrality;
            ThisDataSetSynData(synnum).('Channel_3_CentralitySqrt') = ch3centralitysqrt;
            ThisDataSetSynData(synnum).('Channel_3_WatanabeCentrality') = ch3watanabecentrality;
            
            ThisDataSetSynData(synnum).('Channel_4_SynArea') = ch4synarea;
            ThisDataSetSynData(synnum).('Channel_4_LocNum') = size(thisch4,1);
            ThisDataSetSynData(synnum).('Channel_4_SynDensity') = size(thisch4,1)/ch4synarea;
            ThisDataSetSynData(synnum).('Channel_4_AutoCorr') = ch4AC;
            ThisDataSetSynData(synnum).('Channel_4_Centrality') = ch4centrality;
            ThisDataSetSynData(synnum).('Channel_4_CentralitySqrt') = ch4centralitysqrt;
            ThisDataSetSynData(synnum).('Channel_4_WatanabeCentrality') = ch4watanabecentrality;
            
            ThisDataSetSynData(synnum).('CrossCorr_Ch1_Ch2') = crosscorr12;
            ThisDataSetSynData(synnum).('CrossCorr_Ch1_Ch3') = crosscorr13;
            ThisDataSetSynData(synnum).('CrossCorr_Ch1_Ch4') = crosscorr14;
            ThisDataSetSynData(synnum).('CrossCorr_Ch2_Ch3') = crosscorr23;
            ThisDataSetSynData(synnum).('CrossCorr_Ch2_Ch4') = crosscorr24;
            ThisDataSetSynData(synnum).('CrossCorr_Ch3_Ch4') = crosscorr34;
            
            % append NC params to each synapse in the struct
            fields = fieldnames(NCparams);
            for fieldnum = 1:numel(fields)
                thisfield = fields{fieldnum};
                thisChParams = NCparams.(thisfield);
                thisSynEpsilon = thisChParams(thisChParams(:,1)==synnum,2);
                thisSynMinPts = thisChParams(thisChParams(:,1)==synnum,3);
                ThisDataSetSynData(synnum).([thisfield '_epsilon']) = thisSynEpsilon;
                ThisDataSetSynData(synnum).([thisfield '_minpts']) = thisSynMinPts;
            end
            fprintf('Saved!\nMoving on to nanocluster identification...')
            %%%%% End synaptic analyses
            
            %%%%% Begin nanocluster analyses
            
            %For each channel, we have different params we want to use to
            %define NCs based on dataset. We initialize them per channel here,
            %detect NC and find their peaks for later
            thissyn_withNCtag = [];         %Creating this here for use later
            for channelnum = 1:numel(relevant_channels)
                
                NC_Channel = relevant_channels(channelnum);
                nc_ID = ProteinsInQuestionIDs{NC_Channel};
                protein_to_detect_NCs_in = thissyn(thissyn(:,channelcol)==NC_Channel,:);
                
                % NC parameters from the struct
                nc_epsilon = ThisDataSetSynData(synnum).(['Channel_' num2str(NC_Channel) '_epsilon']);
                nc_minpts = ThisDataSetSynData(synnum).(['Channel_' num2str(NC_Channel) '_minpts']);
                nc_areathresh = outliers.(['Channel_' num2str(NC_Channel)]);
                
                if isempty(protein_to_detect_NCs_in)
                    disp(['Skipping NC detection in ',nc_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset because no locs.']);
                    continue
                end
                
                %Here we actually define NCs.
                NC_of_this_protein =  dbscan(protein_to_detect_NCs_in(:,[xcol ycol]),nc_epsilon,nc_minpts);
                
                %Filters out NCs that we do not want to analyze, either because
                %of low loc number or because they are too big. Often first run we are not
                %filtering out based on area; detect outliers later and
                %rerun
                for ncnum = 1:max(NC_of_this_protein)
                    % Trash NCs that have fewer than 5 locs and move onto
                    % the next one if we do
                    checksizelocs = protein_to_detect_NCs_in(NC_of_this_protein==ncnum,[xcol ycol]);
                    if size(checksizelocs,1) < 5
                        NC_of_this_protein(NC_of_this_protein==ncnum)=-1;
                        continue
                    end
                    
                    % Try to make a shape for the NC; if we can't, then trash the NC
                    % (probably won't happen since that's only an issue if
                    % there are <5 locs, but for safety) and then check its
                    % size and make sure its < areathresh for this channel          
                    thisalphashape = alphaShape(checksizelocs(:,1),checksizelocs(:,2),150/160);
                    [~,ncbound] = boundaryFacets(thisalphashape);
                    if isempty(ncbound)
                       NC_of_this_protein(NC_of_this_protein==ncnum)=-1;
                       continue
                    elseif ~isempty(ncbound)
                        thisncshape = polyshape(ncbound(:,1), ncbound(:,2));
                        areacheck = area(thisncshape);
                        if areacheck > nc_areathresh
                            NC_of_this_protein(NC_of_this_protein==ncnum)=-1;
                        end
                    end

                end
                NC_indtform = [unique(NC_of_this_protein(NC_of_this_protein>0)) [1:(size(unique(NC_of_this_protein(NC_of_this_protein>0)),1))]'];
                for ncind = 1:size(NC_indtform,1)
                    NC_of_this_protein(NC_of_this_protein==NC_indtform(ncind,1)) = NC_indtform(ncind,2);
                end
                
                temp_ncs = [protein_to_detect_NCs_in(:,[xcol ycol]) NC_of_this_protein];
                
                % Determine NC peak so its not calculated twice since it can be
                % variable due to randomization (while the centroid is not)
                ncnum = max(temp_ncs(:,3));
                temparr = [protein_to_detect_NCs_in, NC_of_this_protein, nan(size(NC_of_this_protein,1),2)];
                for thisncnum = 1:ncnum
                    % define peak and add to temparr for later use
                    thisnc = temp_ncs(temp_ncs(:,3)==thisncnum,1:2);
                    thispeak = findDBSCANNCpeak(thisnc,psize,'randomizations',1000);
                    
                    sz = size(thisnc,1);
                    temparr(temparr(:,ncIDcol)==thisncnum,[peakxcol peakycol]) = repmat(thispeak,sz,1); % assign to data
                end
                
                % append temparr to the NCtag variable for use later
                thissyn_withNCtag = [thissyn_withNCtag; temparr];
                
            end
            fprintf('Identified!\n')
            
            %By the end of the above code you have the synapse locs
            %reconstituted with each channel also having the NC designation of
            %each loc appended to it as well as the nc peaks
            
            %This next section starts analyzing things.
            
            % Loop across channels
            for thispiqchannel = 1:length(iq_channels)
                
                % Initialize some things
                proteinInQuestionChannel = iq_channels(thispiqchannel);
                piq_ID = ProteinsInQuestionIDs{proteinInQuestionChannel};
                piqnctotalindex = [];
                
                fprintf('\tAnalyzing %s (%d of %d) nanoclusters...',piq_ID,thispiqchannel,length(iq_channels))
                
                %Grab the locs of the protein in question
                proteinInQuestion = thissyn_withNCtag(thissyn_withNCtag(:,channelcol)==proteinInQuestionChannel,:);
                
                if isempty(proteinInQuestion) % skip if there are no locs
                    disp(['Skipping analysis of ',piq_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset because no locs.']);
                    continue
                end
                
                NC_ProteinInQuestion = proteinInQuestion(:,ncIDcol); %Separate and store NC designations
                n_piqNC = max(NC_ProteinInQuestion);           
                                    
                % Store NC number to the synaptic dataset here in case number is 0 (otherwise skipped and get empty array) 
                if n_piqNC~=-1
                    ThisDataSetSynData(synnum).(['Channel_' num2str(proteinInQuestionChannel) '_NCnum']) = n_piqNC;
                elseif n_piqNC==-1
                    ThisDataSetSynData(synnum).(['Channel_' num2str(proteinInQuestionChannel) '_NCnum']) = 0;
                end

                if n_piqNC == -1 % skip if there are no nanoclusters - we can't do NC-specific analyses if there are
                    disp(['Skipping ',piq_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset because no NCs.']);                     
                    continue
                end
                
                % This next chunk of code randomizes location of NCs
                % BUT - does not adjust the 'hole' left by those NCs, so if you
                % want to do real NC position randomizations, fix that. PD,
                % 11.21.22
                if randflag == 1
                    newpiq = {};
                    newNC_piq = {};
                    newpiqcenters = [];
                    for randrunnum = 1:randruns
                        newpiq{randrunnum} = proteinInQuestion(NC_ProteinInQuestion==-1,1:end-1); %Grabs the nonNC locs
                        newpiq{randrunnum}(:,end+1) = ones(size(newpiq{randrunnum},1),1).*-1; %Gives them a designation of -1
                        for ncnum = 1:max(NC_ProteinInQuestion)
                            newnclocs = [];
                            thisnclocs = proteinInQuestion(NC_ProteinInQuestion==ncnum,1:end-1); %Grab the NC locs for a given NC
                            thisnccenter = [mean(thisnclocs(:,xcol)) mean(thisnclocs(:,ycol))]; %Calculate the center of mass of that NC
                            inbounds = 0;
                            while inbounds == 0 %This While loop randomizes the position of the NC by ensuring its center is within the bounds of the MAGUK
                                ncx = min(thissyn(thissyn(:,channelcol)==3,xcol)) + (max(thissyn(thissyn(:,channelcol)==3,xcol))-min(thissyn(thissyn(:,channelcol)==3,xcol))) .* rand(1,1);
                                ncy = min(thissyn(thissyn(:,channelcol)==3,ycol)) + (max(thissyn(thissyn(:,channelcol)==3,ycol))-min(thissyn(thissyn(:,channelcol)==3,ycol))) .* rand(1,1);
                                synshapetouse = shapeinstructionarray(shapeinstructionarray(:,1)==proteinInQuestionChannel,2);
                                inbounds = inShape(synapseshapestorage{synshapetouse},ncx,ncy);
                                
                                
                                xdiff = ncx-thisnccenter(1);
                                ydiff = ncy-thisnccenter(2);
                                
                                newnclocs(:,1) = thisnclocs(:,xcol)+xdiff;
                                newnclocs(:,2) = thisnclocs(:,ycol)+ydiff;
                                
                            end
                            newpiqcenters(randrunnum,:) = [ncx ncy];
                            thisnclocs(:,2) = newnclocs(:,1);
                            thisnclocs(:,3) = newnclocs(:,2);
                            thisnclocs(:,end+1) = ones(size(thisnclocs,1),1).*ncnum;
                            newpiq{randrunnum} = [newpiq{randrunnum};thisnclocs];
                        end
                        
                        newNC_piq{randrunnum} = newpiq{randrunnum}(:,end);
                    end
                end % end NC position randomization
                
                % get the peaks and centroids out of the data
                ncpiqpeaks = unique(proteinInQuestion(proteinInQuestion(:,ncIDcol)~=-1,[ncIDcol peakxcol peakycol]),'rows');
                ncpiqpeaks = ncpiqpeaks(:,2:3); % crop out the NC number; it was needed to sort correctly by NC # order
                
                % Determine NC peak, shape, and polyshape
                ncpiq = {};
                %ncpiqpeaks = [];
                for piqncnum = 1:max(NC_ProteinInQuestion)
                    thispiqnc = proteinInQuestion(NC_ProteinInQuestion==piqncnum,[xcol ycol]);
                    %ncpiqpeaks(piqncnum,:) = findDBSCANNCpeak(thispiqnc,psize,'randomizations',1000); % calculated above
                    thisncshape = alphaShape(thispiqnc(:,1),thispiqnc(:,2),150/160);
                    [~,ncbound] = boundaryFacets(thisncshape);
                    thisncshape = polyshape(ncbound(:,1), ncbound(:,2));
                    ncpiq{piqncnum} = thisncshape;
                end % end real polyshape calculation
                
                
                % Do the same thing as above but for the random NCs.
                % NOTE THAT as of 5/15/23 ADL the randomized NCs are not going
                % to be precalculated so could ahve different peaks across
                % different parts of this code, since the NC randomization is
                % happening after the peak ID for the real data. Something to
                % fix/adjust later if desired/using randomized data
                if randflag == 1
                    newncpiq = {};
                    newncpiqpeaks = {};
                    for randrunnum = 1:randruns
                        for newpiqncnum = 1:max(newNC_piq{randrunnum})
                            thisnewpiqnc = newpiq{randrunnum}(newNC_piq{randrunnum}==newpiqncnum,[xcol ycol]);
                            newncpiqpeaks{newpiqncnum,1,randrunnum} = findDBSCANNCpeak(thisnewpiqnc,psize,'randomizations',1000);
                            thisncshape = alphaShape(thisnewpiqnc(:,1),thisnewpiqnc(:,2),150/160);
                            [~,ncbound] = boundaryFacets(thisncshape);
                            thisncshape = polyshape(ncbound(:,1), ncbound(:,2));
                            newncpiq{newpiqncnum,1,randrunnum} = thisncshape;
                        end
                    end
                end % end rand polyshape calculation
                
                
                %Here, we use cellfun to calculate the centroid of every
                %polyshape of the real NCs.
                piqnccentroids = {};
                [piqnccentroids(:,1),piqnccentroids(:,2)] = cellfun(@centroid,ncpiq,'UniformOutput',0);
                piqnccentroids = cell2mat(piqnccentroids);
                
                %Here we do the same as above but for the randomly positioned NCs
                if randflag == 1
                    newpiqnccentroids = {};
                    for randrunnum = 1:randruns
                        temp = {};
                        [temp(:,1),temp(:,2)] = cellfun(@centroid,newncpiq(:,:,randrunnum),'UniformOutput',0);
                        newpiqnccentroids{randrunnum} = cell2mat(temp);
                    end
                end % end rand centroid calculation
                
                %Determine center to edge/centrality of nanoclusters relative to reference localizations
                c2eratiorec = [];
                randc2eratiofinal = {};
                % Determine the reference channel based on the shape
                % instruction array
                refchannel = shapeinstructionarray(shapeinstructionarray(:,1)==proteinInQuestionChannel,2);
                clusterlocs = thissyn(thissyn(:,channelcol) == refchannel,[xcol ycol]);
                for piqncnum = 1:n_piqNC
                    
                    thisnccenter = piqnccentroids(piqncnum,:);
                    
                    [centrality, centralitysqrt, watanabecentrality] = Centrality(clusterlocs,thisnccenter,'unitcircle',1);
                    c2eratiorec = [c2eratiorec; centrality, centralitysqrt, watanabecentrality];
                    
                    if randflag == 1
                        randc2eratiorec = []; %Repeat above for random NC
                        for randrunnum = 1:randruns
                            randthisnccenter = newpiqnccentroids{randrunnum}(piqncnum,:);
                            
                            [randcentrality, randcentralitysqrt, randwatanabecentrality] = Centrality(clusterlocs,randthisnccenter,'unitcircle',1);
                            randc2eratiorec = [randc2eratiorec; randcentrality, randcentralitysqrt, randwatanabecentrality];
                            
                        end
                        randc2eratiofinal{piqncnum} = randc2eratiorec;
                    end % end random centrality
                    
                end % end real centrality
                
                
                %At this point we have completed the NC specific analyses for
                %each protein so we store them in a structure for each protein.
                if proteinInQuestionChannel == 1
                    %ThisDataSetSynData(synnum).('Channel_1_NCnum') = n_piqNC;
                    for piqncnum = 1:n_piqNC
                        Channel_1_struct(channel1nctotal).([piq_ID,'_Syn_num']) = synnum;
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_num']) = piqncnum;
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_Area']) = area(ncpiq{piqncnum});
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_Centroid']) = piqnccentroids(piqncnum,:);
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_PeakDensityPoint']) = ncpiqpeaks(piqncnum,:);
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_EllipseCentrality_Real']) = c2eratiorec(piqncnum,1);
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Real']) = c2eratiorec(piqncnum,2);
                        Channel_1_struct(channel1nctotal).([piq_ID,'_NC_WatanabeCentrality_Real']) = c2eratiorec(piqncnum,3);
                        if randflag == 1
                            Channel_1_struct(channel1nctotal).([piq_ID,'_NC_CenterToEdge_Rand']) = randc2eratiofinal{piqncnum};
                            Channel_1_struct(channel1nctotal).([piq_ID,'_NC_EllipseCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,1);
                            Channel_1_struct(channel1nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Rand']) = randc2eratiofinal{piqncnum}(:,2);
                            Channel_1_struct(channel1nctotal).([piq_ID,'_NC_WatanabeCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,3);
                        end
                        piqnctotalindex = [piqnctotalindex; piqncnum channel1nctotal];
                        channel1nctotal = channel1nctotal+1;
                    end
                elseif proteinInQuestionChannel == 2
                    %ThisDataSetSynData(synnum).('Channel_2_NCnum') = n_piqNC;
                    for piqncnum = 1:n_piqNC
                        Channel_2_struct(channel2nctotal).([piq_ID,'_Syn_num']) = synnum;
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_num']) = piqncnum;
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_Area']) = area(ncpiq{piqncnum});
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_Centroid']) = piqnccentroids(piqncnum,:);
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_PeakDensityPoint']) = ncpiqpeaks(piqncnum,:);
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_EllipseCentrality_Real']) = c2eratiorec(piqncnum,1);
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Real']) = c2eratiorec(piqncnum,2);
                        Channel_2_struct(channel2nctotal).([piq_ID,'_NC_WatanabeCentrality_Real']) = c2eratiorec(piqncnum,3);
                        if randflag == 1
                            Channel_2_struct(channel2nctotal).([piq_ID,'_NC_CenterToEdge_Rand']) = randc2eratiofinal{piqncnum};
                            Channel_2_struct(channel2nctotal).([piq_ID,'_NC_EllipseCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,1);
                            Channel_2_struct(channel2nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Rand']) = randc2eratiofinal{piqncnum}(:,2);
                            Channel_2_struct(channel2nctotal).([piq_ID,'_NC_WatanabeCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,3);
                        end
                        piqnctotalindex = [piqnctotalindex; piqncnum channel2nctotal];
                        channel2nctotal = channel2nctotal+1;
                    end
                elseif proteinInQuestionChannel == 3
                    %ThisDataSetSynData(synnum).('Channel_3_NCnum') = n_piqNC;
                    for piqncnum = 1:n_piqNC
                        Channel_3_struct(channel3nctotal).([piq_ID,'_Syn_num']) = synnum;
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_num']) = piqncnum;
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_Area']) = area(ncpiq{piqncnum});
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_Centroid']) = piqnccentroids(piqncnum,:);
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_PeakDensityPoint']) = ncpiqpeaks(piqncnum,:);
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_EllipseCentrality_Real']) = c2eratiorec(piqncnum,1);
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Real']) = c2eratiorec(piqncnum,2);
                        Channel_3_struct(channel3nctotal).([piq_ID,'_NC_WatanabeCentrality_Real']) = c2eratiorec(piqncnum,3);
                        if randflag == 1
                            Channel_3_struct(channel3nctotal).([piq_ID,'_NC_CenterToEdge_Rand']) = randc2eratiofinal{piqncnum};
                            Channel_3_struct(channel3nctotal).([piq_ID,'_NC_EllipseCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,1);
                            Channel_3_struct(channel3nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Rand']) = randc2eratiofinal{piqncnum}(:,2);
                            Channel_3_struct(channel3nctotal).([piq_ID,'_NC_WatanabeCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,3);
                        end
                        piqnctotalindex = [piqnctotalindex; piqncnum channel3nctotal];
                        channel3nctotal = channel3nctotal+1;
                    end
                elseif proteinInQuestionChannel == 4
                    %ThisDataSetSynData(synnum).('Channel_4_NCnum') = n_piqNC;
                    for piqncnum = 1:n_piqNC
                        Channel_4_struct(channel4nctotal).([piq_ID,'_Syn_num']) = synnum;
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_num']) = piqncnum;
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_Area']) = area(ncpiq{piqncnum});
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_Centroid']) = piqnccentroids(piqncnum,:);
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_PeakDensityPoint']) = ncpiqpeaks(piqncnum,:);
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_EllipseCentrality_Real']) = c2eratiorec(piqncnum,1);
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Real']) = c2eratiorec(piqncnum,2);
                        Channel_4_struct(channel4nctotal).([piq_ID,'_NC_WatanabeCentrality_Real']) = c2eratiorec(piqncnum,3);
                        if randflag == 1
                            Channel_4_struct(channel4nctotal).([piq_ID,'_NC_CenterToEdge_Rand']) = randc2eratiofinal{piqncnum};
                            Channel_4_struct(channel4nctotal).([piq_ID,'_NC_EllipseCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,1);
                            Channel_4_struct(channel4nctotal).([piq_ID,'_NC_EllipseCentralitySqrt_Rand']) = randc2eratiofinal{piqncnum}(:,2);
                            Channel_4_struct(channel4nctotal).([piq_ID,'_NC_WatanabeCentrality_Rand']) = randc2eratiofinal{piqncnum}(:,3);
                        end
                        piqnctotalindex = [piqnctotalindex; piqncnum channel4nctotal];
                        channel4nctotal = channel4nctotal+1;
                    end
                end
                fprintf('Saved NC-specific properties!\n\tMoving on to target protein analyses\n')
                %Now that we have stored all the NC properties type analyses,
                %we now conduct analyses comparing the PIQ NCs to TPs
                
                for thistpchannel = 1:length(tp_channels)
                    % initialize a bunch of arrays for outputs
                    AIrec = [];
                    POrec = [];
                    p2prec = [];
                    closesttpnc = {};
                    p2ppairtp = [];
                    ncpiqceoutput = [];
                    ncpiqceoutputcentroid = [];
                    ncpiqisaligned1 = [];
                    ncpiqisaligned1centroid = [];
                    ncpiqisaligned2 = [];
                    ncpiqisaligned2centroid = [];
                    randncpiqoutput = [];
                    newce_realtp_realpiqnc = [];
                    newce_realtp_randpiqNC_median = [];
                    newce_output = [];
                    
                    targetproteinchannel = tp_channels(thistpchannel);
                    tp_ID = TargetProteinIDs{targetproteinchannel};
                    targetProtein = thissyn_withNCtag(thissyn_withNCtag(:,channelcol)==targetproteinchannel,:);
                    fprintf('\t\tAnalyzing %s NC relative to %s (%d of %d)...',piq_ID,tp_ID,thistpchannel,length(tp_channels))
                    
                    if isempty(targetProtein) % skip this loop if there are no locs
                        disp(['Skipping comparing ',piq_ID,' NCs to ',tp_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset because no ',tp_ID,' locs.']);
                        continue
                    end
                    
                    
                    %If there are no NCs in either synapse, we skip it. Because
                    %the AI, PO, p2p all depend on NCs being present in both
                    %protein clusters, we have to put some conditionals
                    %hereafter to account for that. The cross-enrichment only
                    %needs NCs in at least one cluster - this asymmetry is what
                    %complicates the code.
                    NC_TargetProtein = targetProtein(:,ncIDcol);%targetProtein(:,end);
                    
                    if isempty(NC_TargetProtein)
                        n_tpNC = -1;
                    else
                        n_tpNC = max(NC_TargetProtein);
                    end
                    
                    if n_piqNC==-1 && n_tpNC==-1 % skip if no NCs in either synapse as above
                        disp(['Skipping comparing ',piq_ID,' NCs to ',tp_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset because neither have NCs.']);
                        continue
                    end
                    
                    % Get target protein NC peak, shape, polyshape because we
                    % need peaks/centroids for some measures
                    if n_tpNC > - 1
                        
                        % retrieve peaks from previous calculation above to
                        % avoid recalculating given that findDBSCANNCpeak uses
                        % a randomization
                        nctppeaks = unique(targetProtein(targetProtein(:,ncIDcol)~=-1,[ncIDcol peakxcol peakycol]),'rows');
                        nctppeaks = nctppeaks(:,2:3); % crop out the NC number
                        
                        nctp = {};
                        %nctppeaks = []; % NC peaks have been calculated already above, this is leftover
                        for tpncnum = 1:max(NC_TargetProtein)
                            thistpnc = targetProtein(NC_TargetProtein==tpncnum,[xcol ycol]);
                            %nctppeaks(tpncnum,:) = findDBSCANNCpeak(thistpnc,psize,'randomizations',1000);
                            thisncshape = alphaShape(thistpnc(:,1),thistpnc(:,2),150/160);
                            [~,ncbound] = boundaryFacets(thisncshape);
                            thisncshape = polyshape(ncbound(:,1), ncbound(:,2));
                            nctp{tpncnum} = thisncshape;
                        end
                    end % end tp polyshape/peak
                    
                    % get centroids of polyshapes for target protein NC
                    if n_tpNC > -1
                        tpnccentroids = {};
                        [tpnccentroids(:,1),tpnccentroids(:,2)] = cellfun(@centroid,nctp,'UniformOutput',0);
                        tpnccentroids = cell2mat(tpnccentroids);
                    end % end tp centroid
                    
                    % if there are NC in both channels...
                    if n_tpNC>-1 && n_piqNC>-1
                        
                        % get centroid-centroid and p2p distances; as well as
                        % ID of closeset peak and centroid
                        for piqncnum = 1:size(piqnccentroids,1)
                            thispiqnc = piqnccentroids(piqncnum,:);
                            thispiqncpeak = ncpiqpeaks(piqncnum,:);
                            disttotp = pdist2(thispiqnc,tpnccentroids);
                            p2pdist = pdist2(thispiqncpeak,nctppeaks);
                            closesttpnc{piqncnum,1} = thispiqnc; % center of piq NC
                            closesttpnc{piqncnum,2} = tpnccentroids(disttotp==min(disttotp),:) ; %center of closest tp nc
                            closesttpnc{piqncnum,3} = ncpiq{piqncnum}; %piq nc polyshape
                            closesttpnc{piqncnum,4} = nctp{disttotp==min(disttotp)}; %closest tp nc polyshape
                            p2ppairtp(piqncnum,:) = [thispiqncpeak nctppeaks(p2pdist==min(p2pdist),:)]; % xy of piq peak and closet tp peak
                        end % end closest shape/xy calculation
                        if randflag == 1 % repeat on randomized data
                            randclosesttpnc = {};
                            p2ppairtp_rand = {};
                            for randrunnum = 1:randruns
                                for newpiqncnum = 1:size(newpiqnccentroids{randrunnum},1)
                                    thisnewpiqnc = newpiqnccentroids{randrunnum}(newpiqncnum,:);
                                    thispiqncpeak_rand = newncpiqpeaks{newpiqncnum,1,randrunnum};
                                    disttotp = pdist2(thisnewpiqnc,tpnccentroids);
                                    p2pdist_rand = pdist2(thispiqncpeak_rand,nctppeaks);
                                    randclosesttpnc{newpiqncnum,1,randrunnum} = thisnewpiqnc ;
                                    randclosesttpnc{newpiqncnum,2,randrunnum} = tpnccentroids(disttotp==min(disttotp),:);
                                    randclosesttpnc{newpiqncnum,3,randrunnum} = newncpiq{newpiqncnum,1,randrunnum};
                                    randclosesttpnc{newpiqncnum,4,randrunnum} = nctp{disttotp==min(disttotp)};
                                    p2ppairtp_rand{newpiqncnum,1,randrunnum} = [thispiqncpeak_rand nctppeaks(p2pdist_rand==min(p2pdist_rand),:)];
                                end
                            end
                        end % end rand distance/shape closeness calculation
                    end
                    
                    % calculate AI and percent overlap relative to closest tpNC
                    % for real and rand data
                    if n_tpNC>-1 && n_piqNC>-1
                        adjacencyratio = [];
                        pctoverlapout = [];
                        distfact = 5000;
                        for piqncnum = 1:size(closesttpnc,1)
                            
                            if proteinInQuestionChannel ~= targetproteinchannel
                                % calculate adjacency ratio - for the nearest
                                % tp NC, essentially normalize the % overlap to
                                % the size of the clusters
                                v1 = closesttpnc{piqncnum,2}(1)-closesttpnc{piqncnum,1}(1);
                                v2 = closesttpnc{piqncnum,2}(2)-closesttpnc{piqncnum,1}(2);
                                vmag = abs(sqrt((v1^2)+(v2^2)));
                                normv1 = v1/vmag;
                                normv2 = v2/vmag;
                                xt1 = closesttpnc{piqncnum,1}(1)+(distfact*normv1);
                                yt1 = closesttpnc{piqncnum,1}(2)+(distfact*normv2);
                                xt2 = closesttpnc{piqncnum,2}(1)-(distfact*normv1);
                                yt2 = closesttpnc{piqncnum,2}(2)-(distfact*normv2);
                                
                                [in1,~] = intersect(closesttpnc{piqncnum,3},[closesttpnc{piqncnum,1}; xt1 yt1]);
                                [in2,~] = intersect(closesttpnc{piqncnum,4},[closesttpnc{piqncnum,2}; xt2 yt2]);
                                
                                c2c_dist = pdist2([closesttpnc{piqncnum,1}],[closesttpnc{piqncnum,2}]);
                                
                                rad1 = pdist2(in1(1,:), in1(2,:));
                                rad2 = pdist2(in2(1,:), in2(2,:));
                                ratio = c2c_dist/(rad1+rad2);
                            elseif proteinInQuestionChannel == targetproteinchannel
                                ratio = 0;
                            end % end real Adjacency index calculation
                            
                            % calculate % overlap of the piqNC to the nearest
                            % tp NC
                            pctoverlap = area(intersect(closesttpnc{piqncnum,3},closesttpnc{piqncnum,4}))/area(closesttpnc{piqncnum,3});
                            
                            % Repeat for random data
                            if randflag == 1
                                ratiorandrec = [];
                                pctoverlaprandrec = [];
                                for randrunnum = 1:randruns
                                    if proteinInQuestionChannel ~= targetproteinchannel
                                        rv1 = randclosesttpnc{piqncnum,2,randrunnum}(1)-randclosesttpnc{piqncnum,1,randrunnum}(1);
                                        rv2 = randclosesttpnc{piqncnum,2,randrunnum}(2)-randclosesttpnc{piqncnum,1,randrunnum}(2);
                                        rvmag = abs(sqrt((rv1^2)+(rv2^2)));
                                        rnormv1 = rv1/rvmag;
                                        rnormv2 = rv2/rvmag;
                                        xt1rand = randclosesttpnc{piqncnum,1,randrunnum}(1)+(distfact*rnormv1);
                                        yt1rand = randclosesttpnc{piqncnum,1,randrunnum}(2)+(distfact*rnormv2);
                                        xt2rand = randclosesttpnc{piqncnum,2,randrunnum}(1)-(distfact*rnormv1);
                                        yt2rand = randclosesttpnc{piqncnum,2,randrunnum}(2)-(distfact*rnormv2);
                                        
                                        [in1rand,~] = intersect(randclosesttpnc{piqncnum,3,randrunnum},[randclosesttpnc{piqncnum,1,randrunnum}; xt1rand yt1rand]);
                                        [in2rand,~] = intersect(randclosesttpnc{piqncnum,4,randrunnum},[randclosesttpnc{piqncnum,2,randrunnum}; xt2rand yt2rand]);
                                        
                                        c2c_dist_rand = pdist2([randclosesttpnc{piqncnum,1,randrunnum}],[randclosesttpnc{piqncnum,2,randrunnum}]);
                                        rad1rand = pdist2(in1rand(1,:),in1rand(2,:));
                                        rad2rand = pdist2(in2rand(1,:), in2rand(2,:));
                                        thisratiorand = c2c_dist_rand/(rad1+rad2);%(rad1rand+rad2rand); % This is to correct for any oddities of orientation mismatch from real data. Swap with commented values to go back to OG
                                        ratiorandrec = [ratiorandrec;thisratiorand];
                                    elseif proteinInQuestionChannel == targetproteinchannel
                                        ratiorandrec = [ratiorandrec; 0];
                                    end % end rand AI caclulation
                                    pctoverlapthisrand = area(intersect(randclosesttpnc{piqncnum,3},randclosesttpnc{piqncnum,4}))/area(randclosesttpnc{piqncnum,3});
                                    pctoverlaprandrec = [pctoverlaprandrec; pctoverlapthisrand];
                                end % end % overlap calculation
                                
                                avgratiorand = nanmedian(ratiorandrec);
                                adjacencyratio = [adjacencyratio;ratio avgratiorand];
                                
                                avgpctoverlaprand = nanmedian(pctoverlaprandrec);
                                pctoverlapout = [pctoverlapout; pctoverlap avgpctoverlaprand];
                                
                            elseif randflag == 0
                                adjacencyratio = [adjacencyratio;ratio];
                                pctoverlapout = [pctoverlapout; pctoverlap];
                            end % end append real data
                        end % end AI/overlap loop
                        
                        AIrec = [AIrec; adjacencyratio];
                        POrec = [POrec; pctoverlapout];
                        
                        % Determine peak to peak distance for closest NCs
                        thisncrandp2pdist = [];
                        persynp2prec = [];
                        for piqncnum = 1:size(p2ppairtp,1)
                            p2preal = pdist2(p2ppairtp(piqncnum,1:2),p2ppairtp(piqncnum,3:4));
                            if randflag == 1 % repeat for random data
                                for randrunnum = 1:randruns
                                    p2pperrand = pdist2(p2ppairtp_rand{piqncnum,1,randrunnum}(1:2),p2ppairtp_rand{piqncnum,1,randrunnum}(3:4));
                                    thisncrandp2pdist = [thisncrandp2pdist;p2pperrand];
                                end
                                p2prandavg = nanmedian(thisncrandp2pdist);
                                persynp2prec(piqncnum,:) = [p2preal p2prandavg];
                            elseif randflag == 0
                                persynp2prec(piqncnum,:) = p2preal;
                            end
                        end % end p2p calculation
                        p2prec = [p2prec;persynp2prec];
                        
                    end
                    
                    %Now for Enrichment calculations
                    
                    factor = 300;     %density factor to generating the rand cluster
                    proteinInQuestionCE = proteinInQuestion(:,[xcol ycol]);
                    targetProteinCE = targetProtein(:,[xcol ycol]);
                    synshapetouse = shapeinstructionarray(shapeinstructionarray(:,1)==proteinInQuestionChannel,2);
                    rand_cluster_tp = get_cluster_randomized_ADL_PD2D_specifyShape(targetProteinCE, factor, synapseshapestorage{synshapetouse}, 0);
                    
                    for thispiqNC = 1:n_piqNC % loop piq NCs
                        
                        indpiqNC = find(NC_ProteinInQuestion(:,1) == thispiqNC);
                        thisnclocs = proteinInQuestionCE(indpiqNC,:);
                        
                        %~~~~~~~~~~~~~~~CALCULATE DENSITY DISTRIBUTIONS RELATIVE TO NC CENTER~~~~~
                        
                        center = ncpiqpeaks(thispiqNC,:); % get the peak density to use as center
                        distp_piqNCctr = pdist2(targetProteinCE,center);
                        newce_realtp_realpiqnc(thispiqNC,:) = histcounts(distp_piqNCctr, radius)+1;
                        disptprand_piqNCctr = pdist2(rand_cluster_tp,center);
                        disptprand_piqNCctr_hc = histcounts(disptprand_piqNCctr,radius);
                        
                        % This section will make it so that no bin has 0 locs;
                        % randomly redistrubtes locs out of bins so all can
                        % have at least 1.
                        attemptflag = 0;
                        binnum = length(radius)-1;
                        for thisbin = 1:binnum
                            rng('shuffle')
                            loccheck = 0;
                            attemptnum = 1;
                            while loccheck == 0
                                randbin = randi(binnum);
                                loccheck = disptprand_piqNCctr_hc(randbin)>0;
                                attemptnum = attemptnum + 1;
                                if attemptnum > 50
                                    attemptflag = 1;
                                    break
                                end
                            end
                            disptprand_piqNCctr_hc(randbin)=disptprand_piqNCctr_hc(randbin)-1;
                        end
                        if attemptflag == 0
                            disptprand_piqNCctr_hc = disptprand_piqNCctr_hc +1;
                        elseif attemptflag == 1
                            disptprand_piqNCctr_hc = zeros(1,33);
                        end
                        
                        disttp_piqNC = histcounts(distp_piqNCctr, radius) ./ disptprand_piqNCctr_hc .* factor;
                        ncpiqceoutput(thispiqNC,:) = disttp_piqNC;
                        disttp_piqNC(isnan(disttp_piqNC)) = 0;
                        
                        % Determine enrichment index of randomized synapse for
                        % "isAligned" measure baseline
                        EIrandtpdist_piqNC = nan(floor(randruns/2),1);
                        for randrunnum = 1:(floor(randruns/2))
                            synshapetouse = shapeinstructionarray(shapeinstructionarray(:,1)==proteinInQuestionChannel,2);
                            this_rand_cluster_tp = get_cluster_randomized_ADL_PD2D_specifyShape(targetProteinCE, 1, synapseshapestorage{synshapetouse}, 0);
                            thisdisptprand_piqNCctr = pdist2(this_rand_cluster_tp,center);
                            thisdisttp_piqNC = histcounts(thisdisptprand_piqNCctr,radius) ./ disptprand_piqNCctr_hc .* factor;
                            thisdisttp_piqNC(isnan(thisdisttp_piqNC))=0;
                            EIrandtpdist_piqNC(randrunnum,1) = nanmean(thisdisttp_piqNC(2:6));
                            
                        end
                        alignthresholds1 = [(nanmedian(EIrandtpdist_piqNC) + (1.96 * nanstd(EIrandtpdist_piqNC))) (nanmedian(EIrandtpdist_piqNC) - (1.96 * nanstd(EIrandtpdist_piqNC)))];
                        alignthresholds2 = [(nanmean(EIrandtpdist_piqNC) + (1.96 * nanstd(EIrandtpdist_piqNC))) (nanmean(EIrandtpdist_piqNC) - (1.96 * nanstd(EIrandtpdist_piqNC)))];
                        testCE = nanmean(disttp_piqNC(2:6));
                        
                        if testCE > alignthresholds1(1)
                            ncpiqisaligned1(thispiqNC,:) = 1;
                        elseif testCE < alignthresholds1(2)
                            ncpiqisaligned1(thispiqNC,:) = -1;
                        elseif testCE <= alignthresholds1(1) && testCE >= alignthresholds1(2)
                            ncpiqisaligned1(thispiqNC,:) = 0;
                        end
                        if testCE > alignthresholds2(1)
                            ncpiqisaligned2(thispiqNC,:) = 1;
                        elseif testCE < alignthresholds2(2)
                            ncpiqisaligned2(thispiqNC,:) = -1;
                        elseif testCE <= alignthresholds2(1) && testCE >= alignthresholds2(2)
                            ncpiqisaligned2(thispiqNC,:) = 0;
                        end
                        
                        
                        
                        %~~~~~~~~~~~~~~~Repeat relative to centroid~~~~~
                        
                        center2 = piqnccentroids(thispiqNC,:); % use the density peak of the cluster as NC center
                        distp_piqNCctr2 = pdist2(targetProteinCE,center2);
                        newce_realtp_realpiqnc2(thispiqNC,:) = histcounts(distp_piqNCctr2, radius)+1;
                        disptprand_piqNCctr2 = pdist2(rand_cluster_tp,center2);
                        disptprand_piqNCctr_hc2 = histcounts(disptprand_piqNCctr2,radius);
                        
                        % This section will make it so that no bin has 0 locs;
                        % randomly redistrubtes locs out of bins so all can
                        % have at least 1.
                        attemptflag = 0;
                        binnum = length(radius)-1;
                        for thisbin = 1:binnum
                            rng('shuffle')
                            loccheck = 0;
                            attemptnum = 1;
                            while loccheck == 0
                                randbin = randi(binnum);
                                loccheck = disptprand_piqNCctr_hc2(randbin)>0;
                                attemptnum = attemptnum + 1;
                                if attemptnum > 50
                                    attemptflag = 1;
                                    break
                                end
                            end
                            disptprand_piqNCctr_hc2(randbin)=disptprand_piqNCctr_hc2(randbin)-1;
                        end
                        if attemptflag == 0
                            disptprand_piqNCctr_hc2 = disptprand_piqNCctr_hc2 +1;
                        elseif attemptflag == 1
                            disptprand_piqNCctr_hc2 = zeros(1,33);
                        end
                        
                        disttp_piqNC2 = histcounts(distp_piqNCctr2, radius) ./ disptprand_piqNCctr_hc2 .* factor;
                        ncpiqceoutputcentroid(thispiqNC,:) = disttp_piqNC2;
                        disttp_piqNC2(isnan(disttp_piqNC2)) = 0;
                        
                        % Determine enrichment index of randomized synapse for
                        % "isAligned" measure baseline
                        EIrandtpdist_piqNC2 = nan(randruns/2,1);
                        for randrunnum = 1:(randruns/2)
                            synshapetouse = shapeinstructionarray(shapeinstructionarray(:,1)==proteinInQuestionChannel,2);
                            this_rand_cluster_tp2 = get_cluster_randomized_ADL_PD2D_specifyShape(targetProteinCE, 1, synapseshapestorage{synshapetouse}, 0);
                            thisdisptprand_piqNCctr2 = pdist2(this_rand_cluster_tp2,center2);
                            thisdisttp_piqNC2 = histcounts(thisdisptprand_piqNCctr2,radius) ./ disptprand_piqNCctr_hc2 .* factor;
                            thisdisttp_piqNC2(isnan(thisdisttp_piqNC2))=0;
                            EIrandtpdist_piqNC2(randrunnum,1) = nanmean(thisdisttp_piqNC2(2:6));
                            
                        end
                        alignthresholds1c = [(nanmedian(EIrandtpdist_piqNC2) + (1.96 * nanstd(EIrandtpdist_piqNC2))) (nanmedian(EIrandtpdist_piqNC2) - (1.96 * nanstd(EIrandtpdist_piqNC2)))];
                        alignthresholds2c = [(nanmean(EIrandtpdist_piqNC2) + (1.96 * nanstd(EIrandtpdist_piqNC2))) (nanmean(EIrandtpdist_piqNC2) - (1.96 * nanstd(EIrandtpdist_piqNC2)))];
                        testCE = nanmean(disttp_piqNC2(2:6));
                        
                        if testCE > alignthresholds1c(1)
                            ncpiqisaligned1centroid(thispiqNC,:) = 1;
                        elseif testCE < alignthresholds1c(2)
                            ncpiqisaligned1centroid(thispiqNC,:) = -1;
                        elseif testCE <= alignthresholds1c(1) && testCE >= alignthresholds1c(2)
                            ncpiqisaligned1centroid(thispiqNC,:) = 0;
                        end
                        if testCE > alignthresholds2c(1)
                            ncpiqisaligned2centroid(thispiqNC,:) = 1;
                        elseif testCE < alignthresholds2c(2)
                            ncpiqisaligned2centroid(thispiqNC,:) = -1;
                        elseif testCE <= alignthresholds2c(1) && testCE >= alignthresholds2c(2)
                            ncpiqisaligned2centroid(thispiqNC,:) = 0;
                        end
                        
                    end
                    
                    
                    
                    
                    if randflag == 1 % repeat for rand data
                        % note that as of 5/15/23 the rand data will not also
                        % do the centroid based measure, it will just do the
                        % peak.
                        
                        if n_piqNC > 0
                            
                            randpiqcerec = [];
                            
                            for piqncnum = 1:n_piqNC
                                
                                thisNC_realdist_realpos = newce_realtp_realpiqnc(piqncnum,:);
                                newce_realtp_randpiqNC = [];
                                for randrunnum = 1:randruns
                                    randcenter = newncpiqpeaks{piqncnum,:,randrunnum};
                                    distp_randpiqNCctr = pdist2(targetProteinCE,randcenter);
                                    disptprand_randpiqNCctr = pdist2(rand_cluster_tp,randcenter);
                                    disptprand_randpiqNCctr_hc = histcounts(disptprand_piqNCctr,radius);
                                    
                                    % This section will make it so that no bin has 0 locs;
                                    % randomly redistrubtes locs out of bins so all can
                                    % have at least 1.
                                    attemptflag = 0;
                                    binnum = length(radius)-1;
                                    for thisbin = 1:binnum
                                        rng('shuffle')
                                        loccheck = 0;
                                        attemptnum = 1;
                                        while loccheck == 0
                                            randbin = randi(binnum);
                                            loccheck = disptprand_randpiqNCctr_hc(randbin)>0;
                                            attemptnum = attemptnum + 1;
                                            if attemptnum > 50
                                                attemptflag = 1;
                                                break
                                            end
                                        end
                                        disptprand_randpiqNCctr_hc(randbin)=disptprand_randpiqNCctr_hc(randbin)-1;
                                    end
                                    if attemptflag == 0
                                        disptprand_randpiqNCctr_hc = disptprand_randpiqNCctr_hc +1;
                                    elseif attemptflag == 1
                                        disptprand_randpiqNCctr_hc = zeros(1,33);
                                    end
                                    
                                    
                                    randdisttp_piqNC = histcounts(distp_randpiqNCctr, radius) ./ disptprand_randpiqNCctr_hc .* factor;
                                    randpiqcerec(randrunnum,:) = randdisttp_piqNC;
                                    newce_realtp_randpiqNC(randrunnum,:) = histcounts(distp_randpiqNCctr, radius);
                                end
                                
                                newce_realtp_randpiqNC_median(piqncnum,:) = nanmedian(newce_realtp_randpiqNC)+1;
                                newce_output(piqncnum,:) = thisNC_realdist_realpos ./ newce_realtp_randpiqNC_median(piqncnum,:);
                                
                                
                                flip_data = flip(randpiqcerec,2);
                                nanmask=isnan(flip_data);
                                nanmaskwithtrailingvalue = logical([cummin(nanmask,2) ones(size(nanmask,1),1)]);
                                shiftednanmaskwithtrailingvalue=circshift(nanmaskwithtrailingvalue,1,2);
                                actualnanmask = flip(shiftednanmaskwithtrailingvalue(:,1:end-1),2);
                                randpiqcerec(actualnanmask)=9999999999;
                                randpiqcerec(isnan(randpiqcerec))=0;
                                randpiqcerec(isinf(randpiqcerec))=NaN;
                                randpiqcerec(randpiqcerec==9999999999)=NaN;
                                
                                temprandpiqout = nanmedian(randpiqcerec);
                                randncpiqoutput = [randncpiqoutput; temprandpiqout];
                                
                            end
                            
                        elseif n_piqNC == -1
                            temprandpiqout = [];
                            randncpiqoutput = [randncpiqoutput; temprandpiqout];
                        end
                    end
                    
                    % initialize some safeties in case there is no data or an
                    % error
                    if isempty(AIrec)
                        AIrec = NaN(n_piqNC,2);
                    end
                    if isempty(POrec)
                        POrec = NaN(n_piqNC,2);
                    end
                    if isempty(p2prec)
                        p2prec = NaN(n_piqNC,2);
                    end
                    if ~exist('closesttpnc','var') || isempty(closesttpnc)
                        closesttpnc = cell(n_piqNC,4);
                        for ctn = 1:numel(closesttpnc)
                            closesttpnc{ctn} = NaN;
                        end
                    end
                    if ~exist('p2ppairtp','var') || isempty(p2ppairtp)
                        p2ppairtp = NaN(n_piqNC,4);
                    end
                    if isempty(ncpiqceoutput)
                        ncpiqceoutput = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(ncpiqceoutputcentroid)
                        ncpiqceoutputcentroid = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(randncpiqoutput)
                        randncpiqoutput = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(newce_realtp_realpiqnc)
                        newce_realtp_realpiqnc = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(newce_realtp_randpiqNC_median)
                        newce_realtp_randpiqNC_median = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(newce_output)
                        newce_output = NaN(n_piqNC,size(radius,2)-1);
                    end
                    if isempty(targetProteinCE)
                        newce_realtp_realpiqnc = NaN(n_piqNC,size(radius,2)-1);
                        newce_realtp_randpiqNC_median = NaN(n_piqNC,size(radius,2)-1);
                        newce_output = NaN(n_piqNC,size(radius,2)-1);
                    end
                    
                    % store to struct
                    if proteinInQuestionChannel == 1
                        for piqncnum = 1:n_piqNC
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).(['Number_of_',tp_ID,'_NCs_in_this_syn']) = n_tpNC;
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Real']) = AIrec(piqncnum,1);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Real']) = POrec(piqncnum,1);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Real']) = p2prec(piqncnum,1).*psize;
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_center']) = closesttpnc{piqncnum,2};
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_peak']) = p2ppairtp(piqncnum,3:4);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Peak_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutput(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_realpiqnc(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Centroid_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutputcentroid(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1centroid(piqncnum,:);
                            Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2centroid(piqncnum,:);
                            if randflag == 1
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Rand']) = AIrec(piqncnum,2);
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Rand']) = POrec(piqncnum,2);
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Rand']) = p2prec(piqncnum,2).*psize;
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Enrichment_to_',tp_ID,'_Dist']) = randncpiqoutput(piqncnum,:);
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_randpiqNC_median(piqncnum,:);
                                Channel_1_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Coorganization_Quotient_',tp_ID,'_Dist']) = newce_output(piqncnum,:);
                            end
                            
                        end
                    elseif proteinInQuestionChannel == 2
                        for piqncnum = 1:n_piqNC
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).(['Number_of_',tp_ID,'_NCs_in_this_syn']) = n_tpNC;
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Real']) = AIrec(piqncnum,1);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Real']) = POrec(piqncnum,1);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Real']) = p2prec(piqncnum,1).*psize;
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_center']) = closesttpnc{piqncnum,2};
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_peak']) = p2ppairtp(piqncnum,3:4);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Peak_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutput(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_realpiqnc(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Centroid_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutputcentroid(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1centroid(piqncnum,:);
                            Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2centroid(piqncnum,:);
                            if randflag == 1
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Rand']) = AIrec(piqncnum,2);
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Rand']) = POrec(piqncnum,2);
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Rand']) = p2prec(piqncnum,2).*psize;
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Enrichment_to_',tp_ID,'_Dist']) = randncpiqoutput(piqncnum,:);
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_randpiqNC_median(piqncnum,:);
                                Channel_2_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Coorganization_Quotient_',tp_ID,'_Dist']) = newce_output(piqncnum,:);
                            end
                            
                        end
                    elseif proteinInQuestionChannel == 3
                        for piqncnum = 1:n_piqNC
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).(['Number_of_',tp_ID,'_NCs_in_this_syn']) = n_tpNC;
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Real']) = AIrec(piqncnum,1);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Real']) = POrec(piqncnum,1);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Real']) = p2prec(piqncnum,1).*psize;
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_center']) = closesttpnc{piqncnum,2};
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_peak']) = p2ppairtp(piqncnum,3:4);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_Peak_NC_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutput(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_realpiqnc(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Centroid_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutputcentroid(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1centroid(piqncnum,:);
                            Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2centroid(piqncnum,:);
                            if randflag == 1
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Rand']) = AIrec(piqncnum,2);
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Rand']) = POrec(piqncnum,2);
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Rand']) = p2prec(piqncnum,2).*psize;
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Enrichment_to_',tp_ID,'_Dist']) = randncpiqoutput(piqncnum,:);
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_randpiqNC_median(piqncnum,:);
                                Channel_3_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Coorganization_Quotient_',tp_ID,'_Dist']) = newce_output(piqncnum,:);
                            end
                        end
                    elseif proteinInQuestionChannel == 4
                        for piqncnum = 1:n_piqNC
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).(['Number_of_',tp_ID,'_NCs_in_this_syn']) = n_tpNC;
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Real']) = AIrec(piqncnum,1);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Real']) = POrec(piqncnum,1);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Real']) = p2prec(piqncnum,1).*psize;
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_center']) = closesttpnc{piqncnum,2};
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NCs_closest_',tp_ID,'_NC_peak']) = p2ppairtp(piqncnum,3:4);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_Peak_NC_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutput(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Peak_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_realpiqnc(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Real_NC_Centroid_Enrichment_to_',tp_ID,'_Dist']) = ncpiqceoutputcentroid(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Median_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned1centroid(piqncnum,:);
                            Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_NC_Centroid_Mean_isAligned_to_',tp_ID,'_Dist']) = ncpiqisaligned2centroid(piqncnum,:);
                            if randflag == 1
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_AI_to_',tp_ID,'_Rand']) = AIrec(piqncnum,2);
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_PO_to_',tp_ID,'_Rand']) = POrec(piqncnum,2);
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_p2pnm_to_',tp_ID,'_Rand']) = p2prec(piqncnum,2).*psize;
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Enrichment_to_',tp_ID,'_Dist']) = randncpiqoutput(piqncnum,:);
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Rand_NC_Coorg_to_',tp_ID,'_Dist']) = newce_realtp_randpiqNC_median(piqncnum,:);
                                Channel_4_struct(piqnctotalindex(piqnctotalindex(:,1)==piqncnum,2)).([piq_ID,'_Coorganization_Quotient_',tp_ID,'_Dist']) = newce_output(piqncnum,:);
                            end
                        end
                    end
                    
                    fprintf('Saved!\n')
                    %disp(['Done comparing ',piq_ID,' NCs to ',tp_ID,' in Synapse ',num2str(synnum),' of ',ds_ID,' dataset.']);
                    
                    
                end % end tp channel loop
                
                
                
                
                
                
            end % end piq channel loop
            fprintf('\t\tSynapse %d complete!\n',synnum)
            
        catch ME % catch errors and save to the errorhandler struct
            fprintf('\n\t\tOops there was an error on synapse %d, it has been logged in the error handler\n',synnum)
            errorhandler.(['Synapse_' num2str(synnum)]).id = ME.identifier;
            errorhandler.(['Synapse_' num2str(synnum)]).msg = ME.message;
            errorhandler.(['Synapse_' num2str(synnum)]).stk = ME.stack;
        end
        
    end % end synapse loop
    
    
    for piq = 1:length(iq_channels)
        thispiq = iq_channels(piq);
        if thispiq == 1
            piq_ID = ProteinsInQuestionIDs{thispiq};
            ThisDataSetNCData.([piq_ID,'_Data']) = Channel_1_struct;
            disp(['Stored ',piq_ID,' Data']);
        elseif thispiq == 2
            piq_ID = ProteinsInQuestionIDs{thispiq};
            ThisDataSetNCData.([piq_ID,'_Data']) = Channel_2_struct;
            disp(['Stored ',piq_ID,' Data']);
        elseif thispiq == 3
            piq_ID = ProteinsInQuestionIDs{thispiq};
            ThisDataSetNCData.([piq_ID,'_Data']) = Channel_3_struct;
            disp(['Stored ',piq_ID,' Data']);
        elseif thispiq == 4
            piq_ID = ProteinsInQuestionIDs{thispiq};
            ThisDataSetNCData.([piq_ID,'_Data']) = Channel_4_struct;
            disp(['Stored ',piq_ID,' Data']);
        end
    end
    
    AllData.([ds_ID,'_NC_Dataset']) = ThisDataSetNCData;
    AllData.([ds_ID,'_Syn_Dataset']) = ThisDataSetSynData;
    disp(['Stored entire ',ds_ID,' Dataset']);
    
end
disp('All Done!');
toc
warning('on','MATLAB:alphaShape:DupPointsBasicWarnId');

save('FINAL_NC_Analysis.mat','AllData');
if isempty(fieldnames(errorhandler))
    fprintf('There were no errors collected during this analysis!')
else
    fprintf('There were errors collected; see saved errorhandler variable')
    save('ErrorHandler.mat','errorhandler');
end
