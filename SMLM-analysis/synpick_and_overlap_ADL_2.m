%% Isolate potential synapses
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% For all bsn clusters in a region
% 0) First test whether bassoon clusters have all 3 other proteins then...
% 1) Get the associated PSD95 clusters, crop munc and cav to bsn alpha
% shape
% 2) screen putative synapses and flag to redraw boundaries
% 3) save all the keeps to putative_synapses and all the revisits to
% revist_synapses
clear; clc
pathtorequiredfunctions = 'Y:\AAAAAA\Function Repository\Picasso code';
addpath(genpath(pathtorequiredfunctions))
pathtorequiredfunctions = 'Y:\Kaesar exchangePAINT_ADL\Aarons code';
addpath(genpath(pathtorequiredfunctions))

% Grab all week folders
D = dir('Week*');
weekpaths = fullfile({D.folder}',{D.name}');

for thisweek = 1:length(weekpaths) % Loop weeks
    
    S = 'Processing week %d of %d \n';
    fprintf(S,thisweek,length(weekpaths))
    
    % Grab all imaging day folders
    thisweekpath = weekpaths{thisweek};
    imdays = dir(fullfile(thisweekpath,'23*'));
    imdaypaths = fullfile({imdays.folder}',{imdays.name}');
    
    for thisimday = 1:length(imdaypaths) % Loop days
        
        S = 'Processing imaging day %d of %d \n';
        fprintf(S,thisimday,length(imdaypaths))
        
        % Grab all region folders
        thisimdaypath = imdaypaths{thisimday};
        rgns = dir(fullfile(thisimdaypath,'Rgn*'));
        rgnpaths = fullfile({rgns.folder}',{rgns.name}');
        
        for thisrgn = 1:length(rgnpaths) % loop regions
            
            % grab the region paths
            thisrgnpath = rgnpaths{thisrgn};
            analysispath = fullfile(thisrgnpath,'Analysis');
            savepath = fullfile(analysispath,'SynapseProcessing');
            if ~exist(savepath,'dir')
                mkdir(savepath)
            end
            
            S = 'Processing region %d of %d \n';
            fprintf(S,thisrgn,length(rgnpaths))
            
            %%%%%%%%% INITIAL ELIMINATION %%%%%%%%%%%%%%%%%%%%%%
            % load localizations
            [hs,bsn,col] = Omniloader(fullfile(analysispath,'Bsn_link_dbscan_CF_in.hdf5'),'verbose',0);
            [~,psd95,~] = Omniloader(fullfile(analysispath,'PSD95_link_dbscan_CF_in.hdf5'),'verbose',0);
            [~,munc,~] = Omniloader(fullfile(analysispath,'Munc13_link.hdf5'),'verbose',0);
            [~,cav,~] = Omniloader(fullfile(analysispath,'CaV21_link.hdf5'),'verbose',0);
            [~,psdclust,psdclustcol] = Omniloader(fullfile(analysispath,'PSD95_link_dbclusters.hdf5'),'verbose',0);
            
            % test whether bsn clusters have any psd95, munc, and cav
            refloc = bsn;
            testloc = {psd95,munc,cav};
            boxnm = 100;
            psize = 160;
            nlocs = 25;
            [bsnFilter, overlap, nonoverlap] = findPSDoverlap_nlocs(...
                refloc,testloc,col,boxnm,psize,nlocs);
            
            %%%%%%%%%%% DIVIDE INTO KEEP AND REJECT %%%%%%%%%%%%%%%%%%%%%%
            % for each bsn group, associate psd groups with it a la nmdar
            % and build a synapse file, and save rejects to revisit
            bsngroups = unique(bsnFilter(:,col.groups));
            ringwidth = 50; % how much to dilate bsn to search for psd
            
            syns = [];
            ns = 1;
            
            newhs = [hs ',channel,synnum'];
            col2 = getColumns(newhs);
            for thissyn = 1:length(bsngroups)
                
                % grab this bsn group and its psd
                thisbsn = bsnFilter(bsnFilter(:,col.groups) == bsngroups(thissyn),:);
                associatedPSD = associate_groups(thisbsn,psd95,col,psdclust,psdclustcol,ringwidth,boxnm,psize);
                if isempty(associatedPSD)   
                    continue % we don't care about anything with no psd95
                end
                thispsd = psd95(ismember(psd95(:,col.groups),associatedPSD),:);            
                
                % add group if needed, channel, and synnum and concatenate.
                % cav and munc groups are set to 10^6 as garbage
                synarray = [thisbsn repmat(1,size(thisbsn,1),1); ...
                            thispsd repmat(2,size(thispsd,1),1)];
                synnum = repmat(ns,size(synarray,1),1);       
                syns = [syns; synarray synnum];  
                ns = ns + 1;
                       
            end
           
            %%%%%%% SCREEN SYNAPSES FOR A NEED TO REDRAW %%%%%%%%%%%%%%%%%%%%

            revisits = [];
            for thissyn = 1:8:ns-1
                figure('Name', 'Enter n to revisit', 'NumberTitle', 'off','units','normalized','outerposition',[0 0 1 1])
                tl = tiledlayout(2,4);
                
                for g = 0:7
                    
                    if g + thissyn > ns-1 % deal with no multiple of 8 at end
                        continue
                    end
                    % get this bsn and psd
                    thisbsn = syns(syns(:,col2.channel)==1 & syns(:,col2.synnum) == thissyn+g,:);
                    thispsd = syns(syns(:,col2.channel)==2 & syns(:,col2.synnum) == thissyn+g,:);
                    
                    % scatter in appropriate plot
                    nexttile(tl)
                    hold on
                    scatter(thispsd(:,col.x),thispsd(:,col.y),1,'.r')
                    scatter(thisbsn(:,col.x),thisbsn(:,col.y),1,'.b')
                    axis equal
                    title(num2str(thissyn+g))
                        
                end
                hl = legend({'psd95','bsn'});
                hl.Layout.Tile = 'East';
                revisit = inputdlg('Enter numbers to revisit'); % request to keep
                % generate revisit list
                if ~isempty(revisit{1})
                    t = strsplit(revisit{1},' ')';
                    revisits = [revisits; str2double(t)];
                end
                
                close
            end
                         
            %%%%%%%%%%%%% SAVE FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            %%%%%%%%% SAVE PUTATIVE GOOD SYNAPSES
            RPTP(syns,newhs,fullfile(analysispath,'CaV21_link'),'null',...
                'changename','potential_synapses','changepath',savepath)
            
            save(fullfile(savepath,'SynsToRevisit'),'revisits')
            
        end %rgn
        
        fprintf('Finished imaging day %d \n',thisimday)
        
    end %imday
    
    fprintf('Finished week %d \n', thisweek)
    
end %week

fprintf('Processing is finished!')

%% Redraw borders on revisit synapses
% grab the revisit file from each region, and then manually redraw borders
% of bsn and psd clusters based on psd/bsn image. Then check as above whether the
% redrawn bsn/psd meets the aspect ratio and % overlap criteria and if so
% keep it. Finally reload the potential synapses file and append the
% reclaimed synapses to it and save as _for_QC. Note that psd95 kept here
% is the psd95 that is drawn around

clear; clc;
% Grab all week folders
D = dir('Week*');
weekpaths = fullfile({D.folder}',{D.name}');

for thisweek = 1:length(weekpaths) % Loop weeks
    
    S = 'Processing week %d of %d \n';
    fprintf(S,thisweek,length(weekpaths))
    
    % Grab all imaging day folders
    thisweekpath = weekpaths{thisweek};
    imdays = dir(fullfile(thisweekpath,'23*'));
    imdaypaths = fullfile({imdays.folder}',{imdays.name}');
    
    for thisimday = 1:length(imdaypaths) % Loop days
        
        S = 'Processing imaging day %d of %d \n';
        fprintf(S,thisimday,length(imdaypaths))
        
        % Grab all region folders
        thisimdaypath = imdaypaths{thisimday};
        rgns = dir(fullfile(thisimdaypath,'Rgn*'));
        rgnpaths = fullfile({rgns.folder}',{rgns.name}');
        
        for thisrgn = 1:length(rgnpaths) % loop regions
            
            % grab the region paths
            thisrgnpath = rgnpaths{thisrgn};
            analysispath = fullfile(thisrgnpath,'Analysis');
            savepath = fullfile(analysispath,'SynapseProcessing');
            if ~exist(savepath,'dir')
                error(['You dont have a path to the synapse processing folder for ' thisrgnpath])
            end
            
            S = 'Processing region %d of %d \n';
            fprintf(S,thisrgn,length(rgnpaths))
            
            % load localizations
            if exist(fullfile(savepath,'potential_synapses.hdf5'),'file')
                [hs,loc,col] = Omniloader(fullfile(savepath,'potential_synapses.hdf5'),'verbose',0);
                [~,munc,~] = Omniloader(fullfile(analysispath,'Munc13_link.hdf5'),'verbose',0);
                [~,cav,~] = Omniloader(fullfile(analysispath,'CaV21_link.hdf5'),'verbose',0);
                load(fullfile(savepath,'SynsToRevisit.mat'));
            else
                continue
            end
            
            nsyns = max(loc(:,col.synnum));
            boxnm = 100;
            psize = 160;
            ns = 1;
            syns = [];
          
            warning('off','MATLAB:polyshape:repairedBySimplify') % suppress a polyshape warning in getPSDoverlap
            for i = 1:nsyns
                
                thissyn = loc(loc(:,col.synnum)==i,:);
                
                if ~ismember(i,revisits)
                    
                    thisbsn = thissyn(thissyn(:,col.channel)==1,:);
                    thispsd = thissyn(thissyn(:,col.channel)==2,:);
                
                elseif ismember(i,revisits) 
                
                    figure
                    hold on
                    axis equal
                    scatter(thissyn(thissyn(:,col.channel)==1,col.x), thissyn(thissyn(:,col.channel)==1,col.y), '.r');
                    scatter(thissyn(thissyn(:,col.channel)==2,col.x), thissyn(thissyn(:,col.channel)==2,col.y), '.b');
                    legend({'psd95','bsn'})
                    title([num2str(i),' of ',num2str(nsyns)]);
                    answer = input('0 to trash, 1 to redraw');

                    if answer == 0
                        close all
                        continue

                    elseif answer == 1
                        % draw new borders and crop psd and bsn to that space
                        h = drawfreehand;
                        syntf = inROI(h,thissyn(:,col.x), thissyn(:,col.y));
                        croppedsyn = thissyn(syntf,:);
                        thisbsn = croppedsyn(croppedsyn(:,col.channel)==1,:);
                        thispsd = croppedsyn(croppedsyn(:,col.channel)==2,:);
                        close all
                        
                    end
                    
                else
                    error('You had something that wasnt a keep or revisit, your syn variable is screwed up')
                end
          
                % calculate aspect ratio and area overlap
                [ratio, ~, ~] = getPSDaxes_alphaRad_ADL(thisbsn(:,[col.x col.y]),'alphaRadius',2);
                area1 = getPSDoverlap(thisbsn(:,[col.x col.y]), thispsd(:,[col.x col.y]));
                area2 = getPSDoverlap(thispsd(:,[col.x col.y]), thisbsn(:,[col.x col.y]));
                maxarea = max(area1,area2);
                    
                % crop data and concatenate synapse files
                if ratio <= 2.0 && maxarea >= 0.7 % for good synapses
                    
                    % crop the cav and munc to inside the bsn shape
                    bsnshp = alphaShape(thisbsn(:,col.x), thisbsn(:,col.y), 'HoleThreshold', 100000);
                    if area(bsnshp) < pi*(20/psize)^2 || area(bsnshp) > pi*(500/psize)^2
                        continue % eliminate anything ~ area of 1 munc13 NC (radius 15 nm-ish) or > 500 nm radius (1000 nm diameter)
                    end
                    box_minx = min(thisbsn(:,col.x)) - boxnm/psize;
                    box_maxx = max(thisbsn(:,col.x)) + boxnm/psize;
                    box_miny = min(thisbsn(:,col.y)) - boxnm/psize;
                    box_maxy = max(thisbsn(:,col.y)) + boxnm/psize;
                    inrectmunc = munc(munc(:,col.x) <= box_maxx & munc(:,col.x) >= box_minx & ...
                        munc(:,col.y) <= box_maxy & munc(:,col.y) >= box_miny,:);
                    inrectcav =  cav(cav(:,col.x) <= box_maxx & cav(:,col.x) >= box_minx & ...
                        cav(:,col.y) <= box_maxy & cav(:,col.y) >= box_miny,:);
                    thismunc = inrectmunc(inShape(bsnshp,inrectmunc(:,[col.x col.y])),:);
                    thiscav = inrectcav(inShape(bsnshp,inrectcav(:,[col.x col.y])),:);
                                   
                    % add group if needed, channel, and synnum and concatenate
                    synarray = [thisbsn(:,1:end-2) repmat(1,size(thisbsn,1),1); ...
                                thispsd(:,1:end-2) repmat(2,size(thispsd,1),1); ...
                                thismunc repmat(10^6,size(thismunc,1),1) repmat(3,size(thismunc,1),1); ...
                                thiscav repmat(10^6,size(thiscav,1),1) repmat(4,size(thiscav,1),1)];
                    synnum = repmat(ns,size(synarray,1),1);
                    syns = [syns; synarray synnum];
                    ns = ns + 1;
                else % if it doesn't pass skip it
                    continue
                end
            end   
            warning('on','MATLAB:polyshape:repairedBySimplify')      

            %%%%%%% SAVE THE FINAL SYNAPSE FILE
            RPTP(syns,hs,fullfile(savepath,'potential_synapses'),'null',...
                'changename','0_all_synapses_for_QC','changepath',savepath)
            % save individual channels
            RPTP(syns(syns(:,col.channel)==1,:),...
                hs,fullfile(analysispath,'Bsn_link_dbscan_CF_in'),'null',...
                'changename','0_Bsn_synapses_for_QC','changepath',savepath)
            RPTP(syns(syns(:,col.channel)==2,:),...
                hs,fullfile(analysispath,'PSD95_link_dbscan_CF_in'),'null',...
                'changename','0_PSD95_synapses_for_QC','changepath',savepath)
            RPTP(syns(syns(:,col.channel)==3,:),...
                hs,fullfile(analysispath,'Munc13_link'),'null',...
                'changename','0_Munc13_synapses_for_QC','changepath',savepath)
            RPTP(syns(syns(:,col.channel)==4,:),...
                hs,fullfile(analysispath,'CaV21_link'),'null',...
                'changename','0_CaV21_synapses_for_QC','changepath',savepath)
            
        end %rgn
        
        fprintf('Finished imaging day %d \n',thisimday)
        
    end %imday
    
    fprintf('Finished week %d \n', thisweek)
    
end %week

fprintf('Processing is finished!')
%% % plot a random synapse if needed.
figure
synnum = 1;
thisbsn = syns(syns(:,col.channel)==1 & syns(:,col.synnum)==synnum,:);
thispsd = syns(syns(:,col.channel)==2 & syns(:,col.synnum)==synnum,:);
hold on
scatter(thispsd(:,col.x),thispsd(:,col.y),'.r')
scatter(thisbsn(:,col.x),thisbsn(:,col.y),'.b')
axis equal


%% Review potential synapses CaV and Munc and reject if obviously bad
% plot the bsn alpha shape along with cav and munc locs. enter numbers to
% reject, then save hdf5 and synapse.txt files to appropriate locations as 
% final synapses that will go into analysis.

clear; clc;
% Grab all week folders
D = dir('Week*');
weekpaths = fullfile({D.folder}',{D.name}');

syncount = 1;
synpath = fullfile(pwd,'Synapses');
for thisweek = 1:length(weekpaths) % Loop weeks
    
    S = 'Processing week %d of %d \n';
    fprintf(S,thisweek,length(weekpaths))
    
    % Grab all imaging day folders
    thisweekpath = weekpaths{thisweek};
    imdays = dir(fullfile(thisweekpath,'23*'));
    imdaypaths = fullfile({imdays.folder}',{imdays.name}');
    
    for thisimday = 1:length(imdaypaths) % Loop days
        
        S = 'Processing imaging day %d of %d \n';
        fprintf(S,thisimday,length(imdaypaths))
        
        % Grab all region folders
        thisimdaypath = imdaypaths{thisimday};
        rgns = dir(fullfile(thisimdaypath,'Rgn*'));
        rgnpaths = fullfile({rgns.folder}',{rgns.name}');
        
        for thisrgn = 1:length(rgnpaths) % loop regions
            
            % grab the region paths
            thisrgnpath = rgnpaths{thisrgn};
            analysispath = fullfile(thisrgnpath,'Analysis');
            savepath = fullfile(analysispath,'SynapseProcessing');
            if ~exist(savepath,'dir')
                error(['You are missing a savepath in ' thisrgnpath])
            end
            
            S = 'Processing region %d of %d \n';
            fprintf(S,thisrgn,length(rgnpaths))
            
            % load localizations
            if exist(fullfile(savepath,'0_all_synapses_for_QC.hdf5'),'file')
                [hs,loc,col] = Omniloader(fullfile(savepath,'0_all_synapses_for_QC.hdf5'),'verbose',0);           
            else
                continue
            end
            
            synn = max(loc(:,col.synnum));
            garbagecollector = [];
            for thissyn = 1:6:synn
                figure('Name', 'Enter n to trash', 'NumberTitle', 'off','units','normalized','outerposition',[0 0 1 1])
                tl = tiledlayout(3,2);
                
                for g = 0:5
                    
                    if g + thissyn > synn % deal with no multiple of 8 at end
                        continue
                    end
                    % get this bsn and psd
                    thisbsn = loc(loc(:,col.channel)==1 & loc(:,col.synnum)==thissyn+g,:);
                    thispsd = loc(loc(:,col.channel)==2 & loc(:,col.synnum)==thissyn+g,:);
                    thismunc = loc(loc(:,col.channel)==3 & loc(:,col.synnum)==thissyn+g,:);
                    thiscav = loc(loc(:,col.channel)==4 & loc(:,col.synnum)==thissyn+g,:); 
                    
                    % scatter in appropriate plot
                    t2 = tiledlayout(tl,1,2,'Padding','tight','TileSpacing','none');
                    t2.Layout.Tile = g+1;
                    t2.Title.String = ['Synapse ' num2str(thissyn+g)];
                    t2.Title.FontWeight = 'bold';
                    ax1 = nexttile(t2);
                    hold on
                    scatter(thismunc(:,col.x),thismunc(:,col.y),10,'om','filled')
                    scatter(thiscav(:,col.x),thiscav(:,col.y),10,'og','filled')
                    daspect([1 1 1])
                    pbaspect([1 1 1])
                    title('Munc (m), Cav (g)')
                    
                    ax2 = nexttile(t2);
                    hold on
                    scatter(thispsd(:,col.x),thispsd(:,col.y),10,'or','filled')
                    scatter(thisbsn(:,col.x),thisbsn(:,col.y),10,'ob','filled')
                    daspect([1 1 1])
                    pbaspect([1 1 1])
                    title('PSD (r), Bsn (b)')
                    
                    linkaxes([ax1 ax2],'xy')
                        
                end

                garbage = inputdlg('Enter numbers to reject'); % request to keep
                % generate garbage list
                if ~isempty(garbage{1})
                    t = strsplit(garbage{1},' ')';
                    garbagecollector = [garbagecollector; str2double(t)];
                end
                
                close
            end
            
            keepsyns = loc(~ismember(loc(:,col.synnum),garbagecollector),:); 
            
            % also renumber the kept synapses since we removed some
            uniquekeepsyns = unique(keepsyns(:,col.synnum));
            tempkeepsyns = [];
            for r = 1:length(uniquekeepsyns) 
                tempsyn = keepsyns(keepsyns(:,col.synnum) == uniquekeepsyns(r),:);
                tempsyn(:,col.synnum) = r;
                tempkeepsyns = [tempkeepsyns; tempsyn];  
            end 
            keepsyns = tempkeepsyns;
           
            if ~isempty(keepsyns)
                %%%%%%% SAVE THE FINAL SYNAPSE FILE
                % save all synapses
                RPTP(keepsyns,hs,fullfile(savepath,'0_all_synapses_for_QC'),'null',...
                    'changename','1_all_synapses_for_analysis')
                save(fullfile(savepath,'trashedSynapses.mat'),'garbagecollector')

                % save individual channels
                RPTP(keepsyns(keepsyns(:,col.channel)==1,:),...
                    hs,fullfile(savepath,'0_Bsn_synapses_for_QC'),'null','changename','1_Bsn_synapses_for_analysis')
                RPTP(keepsyns(keepsyns(:,col.channel)==2,:),...
                    hs,fullfile(savepath,'0_PSD95_synapses_for_QC'),'null','changename','1_PSD95_synapses_for_analysis')
                RPTP(keepsyns(keepsyns(:,col.channel)==3,:),...
                    hs,fullfile(savepath,'0_Munc13_synapses_for_QC'),'null','changename','1_Munc13_synapses_for_analysis')
                RPTP(keepsyns(keepsyns(:,col.channel)==4,:),...
                    hs,fullfile(savepath,'0_CaV21_synapses_for_QC'),'null','changename','1_CaV21_synapses_for_analysis')                     

                % plot and save a figure
                figure
                hold on
                scatter(keepsyns(keepsyns(:,col.channel)==2,col.x),keepsyns(keepsyns(:,col.channel)==2,col.y),1,'.b');
                scatter(keepsyns(keepsyns(:,col.channel)==1,col.x),keepsyns(keepsyns(:,col.channel)==1,col.y),1,'.g');
                scatter(keepsyns(keepsyns(:,col.channel)==3,col.x),keepsyns(keepsyns(:,col.channel)==3,col.y),1,'.c');
                scatter(keepsyns(keepsyns(:,col.channel)==4,col.x),keepsyns(keepsyns(:,col.channel)==4,col.y),1,'.m')
                axis equal
                legend({'PSD95','Bsn','Munc13','CaV21'})
                savefig(gcf,fullfile(savepath,'2AnalyzedSynapses.fig'))
                close all

                % also save synapse.txt file to parent directory
                bing = size(keepsyns,1);
                keepsyns2 = [repmat(thisweek,bing,1) repmat(thisimday,bing,1) repmat(thisrgn,bing,1) keepsyns];
                writematrix(keepsyns2,fullfile(synpath,['synapse' num2str(syncount) '.txt']),'delimiter','tab');
                syncount = syncount+1;
            end
            
                      
        end %rgn
        
        fprintf('Finished imaging day %d \n',thisimday)
        
    end %imday
    
    fprintf('Finished week %d \n', thisweek)
    
end %week

fprintf('Processing is finished!')
