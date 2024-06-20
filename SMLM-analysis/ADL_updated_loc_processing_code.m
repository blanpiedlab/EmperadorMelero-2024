%% Process Kaesar 4-color Exchange-PAINT images through all localization and filtering steps
% Folder structure
% >Cell week
%    >Imaging day
%        >bead images
%        >Rgn 1 - n
%           >with nd2 files, or with existing _left/right_locs files
% Requires picassosr 0.6.0 installed in python and linked to windows cmd via environment variable
% Must update __main__.py in the picassosr package with __main__0pt6.py here.
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine

clear; clc;
% Add required functions
pathtorequiredfunctions = 'Y:\AAAAAA\Function Repository\Picasso code';
addpath(genpath(pathtorequiredfunctions))

%%%%%%%%OVERALL PARAMS
pixelsize = 160;

%%%%%%%%DUVCALI PARAMS
% Localization
boxsize = 7;
gradient = 20000;
fitmethod = 'lq-gpu'; 
drift = 0; % set to 0 to do independent RCC below
camera = 'ixon';
% DUV cali
direction = 'l2r';
splitpoint = 256;

%%%%%%%%LOC PARAMS
% Localization
reloc = 0; % set to 1 to also localize if the nd2 file exists
boxsize_l = 7;
leftgradient = 9000;
rightgradient = 9000;
leftroi = '0 0 300 256'; %ymin xmin ymax ymax
rightroi = '0 256 300 512';
fitmethod_l = 'lq-gpu'; 
drift_l = 0; % set to 0 to do independent RCC below
camera_l = 'ixon';
%Drift correct
seg = 1000;

%%%%%%%FILTERING PARAMS
% Link
radius = 0.3;
darkframes = 5; 
% Filtering
maxerror = 20/pixelsize;
sigmahi = 1.6; 
sigmalo = 0.3;
% DBSCAN
eps = 48/pixelsize;
minpts = 10;
%ClusterFilter
std_range = [2500 11000];



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
        
        %%%%%%%% SET UP FOR DUV CALI AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        califolder = fullfile(thisimdaypath,'DUVcali_files');
        if ~exist(califolder,'dir')
            
            S = 'Calculating DUV calibration \n';
            fprintf(S)
            
            % grab the bead files and move to a subfolder
            califiles = dir(fullfile(thisimdaypath,'*beads*.nd2'));
            calipaths = fullfile({califiles.folder}',{califiles.name}');
            mkdir(califolder);
            for i = 1:size(calipaths,1)
                [~] = movefile(calipaths{i},califolder);
            end
            califiles = dir(fullfile(califolder,'*beads*.nd2'));
            calipaths = fullfile({califiles.folder}',{califiles.name}');

            % localize DUVcali files
            for i = 1:size(calipaths,1)
                impath = calipaths{i};
                picassoLocalizev060(impath,boxsize,fitmethod,gradient,drift,camera,'verbose',0);
            end

            % join all cali loc files
            calilocfiles = dir(fullfile(califolder,'*locs.hdf5'));
            calilocpaths = fullfile({calilocfiles.folder}',{calilocfiles.name}');
            calijoinpath = picassoJoin(calilocpaths);

            %%%%%%%%%%%%%%% GENERATE TFORM FOR DUV CORRECTION %%%%%%%%%%%%
            ADLmakeTform_poly_bidirectional(calijoinpath, pixelsize, direction, ...
                'Splitpoint', splitpoint,'saveinplace',1);
            savefig(gcf,fullfile(califolder,'duvcali.fig'))
            close all

            % move cali file to region subfolders
            tformpath = fullfile(califolder,['tform' direction '.mat']);
            cellfun(@(x) copyfile(tformpath,x),rgnpaths,'uni',0);
            
        end % cali
        

        %%%%%%%%%%%%%%% PROCESS REGIONS %%%%%%%%%%%%%%%%%%%%%%%%
        
        for thisrgn = 1:length(rgnpaths) % loop regions
            
            % grab the region paths
            thisrgnpath = rgnpaths{thisrgn};
            
            S = 'Processing region %d of %d \n';
            fprintf(S,thisrgn,length(rgnpaths))
            
            % check for tform
            if ~exist(fullfile(thisrgnpath,['tform' direction '.mat']),'file')
                try
                    tformpath = fullfile(califolder,['tform' direction '.mat']);
                    copyfile(tformpath,thisrgnpath)
                catch
                    error(['You have tform issues in ' thisrgnpath ', figure it out and try again'])
                end
            end
            
            %%%%%%%%% LOCALIZE, RE-JOIN, UNDRIFT, DUV CORRECT EACH ROUND %%%%%%%%% 
            for round = 1:2 % loop exchange rounds
                
                S = 'Localize, join, drift and DUV correct round %d \n';
                fprintf(S,round)
                
                % grab the image path
                impath = fullfile(thisrgnpath,['Round' num2str(round) '_.nd2']);

                % localize the left and right side using gradients defined
                % in settings, or else join existing loc files and proceed
                jointemp = cell(2,1);
                if exist(impath,'file') && reloc == 1  %localize left and right of im
                    [jointemp{1}, ~] = picassoLocalizev060(impath,boxsize_l,fitmethod_l,...
                        leftgradient,drift_1,camera_l,'roi',leftroi,'suffix','_left','verbose',0);
                    [jointemp{2}, ~] = picassoLocalizev060(impath,boxsize_l,fitmethod_l,...
                        rightgradient,drift_l,camera_l,'roi',rightroi,'suffix','_right','verbose',0);
                elseif reloc == 0 % find and use existing loc files
                    jointemp{1,1} = fullfile(thisrgnpath,['Round' num2str(round) '_left_locs.hdf5']);
                    jointemp{2,1} = fullfile(thisrgnpath,['Round' num2str(round) '_right_locs.hdf5']);
                    if ~exist(jointemp{1},'file') || ~exist(jointemp{2},'file')
                        error(['You are missing a loc file in ' thisrgnpath])
                    end
                elseif ~exist(impath,'file') && reloc == 1
                    error(['You are missing an nd2 file in ' thisrgnpath])
                end
                
                % join both sides to one file
                joinlocpath = picassoJoin(jointemp,'keep',1);

                % undrift
                [driftpath,~,~] = picassoUndrift(joinlocpath,'seg',seg,'verbose',0);

                % DUV correct
                [hs,totallocs] = ADLcorrectDUV_poly_bidirectional(driftpath,pixelsize,...
                    'SplitPoint',splitpoint,'dontsave',1,'figflag',1,'loadinplace',1);
                col = getColumns(hs);
                savefig(gcf,fullfile(thisrgnpath,['DUV_round_' num2str(round) '.fig']))
                close all
                
                % separate back to left and right and re zero if necessary
                colstokeep = 1:size(totallocs,2)~=col.channel;
                locleft = totallocs(totallocs(:,col.channel)==1,colstokeep);
                locright = totallocs(totallocs(:,col.channel)==2,colstokeep);

                if strcmp(direction,'l2r')==1
                    locleft(:,col.x) = locleft(:,col.x)-splitpoint;
                    locright(:,col.x) = locright(:,col.x)-splitpoint;
                end

                % Resave separate channels
                hs = hs(1:strfind(hs,',Channel')-1);           
                [thisfolder,~,~] = fileparts(driftpath);
                oldname = erase(driftpath,'.hdf5');
                leftname = ['Round' num2str(round) '_left_locs_undrift_DUV'];
                rightname = ['Round' num2str(round) '_right_locs_undrift_DUV'];
                RPTP(locleft, hs, oldname, 'null', 'changename', leftname);
                RPTP(locright, hs, oldname, 'null', 'changename', rightname);  

        
            end % imaging round
        
            %%%%%%%% XCORR EACH IMAGE TO REFERENCE %%%%%%%%%%%%%%%%%%
            S = 'Cross correlating...\n';
            fprintf(S)
            
            fileR1L = fullfile(thisrgnpath,'Round1_left_locs_undrift_DUV.hdf5');
            fileR1R = fullfile(thisrgnpath,'Round1_right_locs_undrift_DUV.hdf5');
            fileR2L = fullfile(thisrgnpath,'Round2_left_locs_undrift_DUV.hdf5');
            fileR2R = fullfile(thisrgnpath,'Round2_right_locs_undrift_DUV.hdf5');
            [hs,R1L,col] = Omniloader(fileR1L,'verbose',0);
            [~,R1R,~] = Omniloader(fileR1R,'verbose',0);        
            [~,R2L,~] = Omniloader(fileR2L,'verbose',0);
            [~,R2R,~] = Omniloader(fileR2R,'verbose',0);
                  
            shiftV = PAINTshift(R1R(:,[col.x col.y]), R1L(:,[col.x col.y]), 5, pixelsize, 50, [], [0 35]);
            R1L(:,col.x) = R1L(:,col.x)+shiftV(1);
            R1L(:,col.y) = R1L(:,col.y)+shiftV(2);
            
            shiftV = PAINTshift(R1R(:,[col.x col.y]), R2L(:,[col.x col.y]), 5, pixelsize, 50, [], [0 35]);
            R2L(:,col.x) = R2L(:,col.x)+shiftV(1);
            R2L(:,col.y) = R2L(:,col.y)+shiftV(2);
            
            shiftV = PAINTshift(R1R(:,[col.x col.y]), R2R(:,[col.x col.y]), 5, pixelsize, 50, [], [0 35]);
            R2R(:,col.x) = R2R(:,col.x)+shiftV(1);
            R2R(:,col.y) = R2R(:,col.y)+shiftV(2);
            
            RPTP(R1L,hs,erase(fileR1L,'.hdf5'),'_xcorr')
            RPTP(R1R,hs,erase(fileR1R,'.hdf5'),'_xcorr')
            RPTP(R2L,hs,erase(fileR2L,'.hdf5'),'_xcorr')
            RPTP(R2R,hs,erase(fileR2R,'.hdf5'),'_xcorr')
            
            % plot and save a figure
            figure
            hold on
            scatter(R1L(:,col.x),R1L(:,col.y),1,'.b')
            scatter(R2L(:,col.x),R2L(:,col.y),1,'.g');
            scatter(R1R(:,col.x),R1R(:,col.y),1,'.c');
            scatter(R2R(:,col.x),R2R(:,col.y),1,'.m');
            legend({'psd95','bsn','munc','cav'})
            axis equal
            savefig(gcf,fullfile(thisrgnpath,'CrossCorr_corrected.fig'));
            close all
            
            %%%%%%%% FILTER AND LINK, MOVE FILES TO ANALYSIS FOLDER %%%%%%
            % create analysis folder
            analysispath = fullfile(thisrgnpath,'Analysis');
            if ~exist(analysispath,'dir')
                mkdir(analysispath)
            end  
            locs = {R1L;R1R;R2L;R2R};
            locpaths = {fileR1L;fileR1R;fileR2L;fileR2R};
            
            for thisloc = 1:length(locs) % loop proteins
                
                S = 'Filtering and linking file %d of %d \n';
                fprintf(S,thisloc,length(locs))
                
                % grab these locs
                loc = locs{thisloc};
                thislocpath = locpaths{thisloc};
                
                % filter on sigma
                loc = loc(loc(:,col.sigmax)<sigmahi & loc(:,col.sigmax)>sigmalo & loc(:,col.sigmay)<sigmahi & loc(:,col.sigmay)>sigmalo,:);

                % filter on mode of photon distribution
                ph = histogram(loc(:,col.phot),'BinWidth',100,'BinCountsMode','auto');
                photonmode = ph.BinEdges(ph.BinCounts == max(ph.BinCounts));
                close all
                loc = loc(loc(:,col.phot)>photonmode,:);

                % filter on error
                loc = loc(loc(:,col.lpx)<maxerror & loc(:,col.lpy)<maxerror,:);
                
                % save filter
                RPTP(loc,hs,erase(thislocpath,'.hdf5'),'_xcorr_filter')
          
                % link 
                linkpath = strrep(thislocpath,'.hdf5','_xcorr_filter.hdf5');
                picassoLink(linkpath, radius, darkframes);
                
                % move and rename files to analysis folder
                movepath = strrep(linkpath,'.hdf5','_link');
                if contains(movepath,'Round1_left')
                    copyfile([movepath '.hdf5'],fullfile(analysispath,'PSD95_link.hdf5'));
                    copyfile([movepath '.yaml'],fullfile(analysispath,'PSD95_link.yaml'));
                elseif contains(movepath,'Round1_right')
                    copyfile([movepath '.hdf5'],fullfile(analysispath,'Munc13_link.hdf5'));
                    copyfile([movepath '.yaml'],fullfile(analysispath,'Munc13_link.yaml'));
                elseif contains(movepath,'Round2_left')
                    copyfile([movepath '.hdf5'],fullfile(analysispath,'Bsn_link.hdf5'));
                    copyfile([movepath '.yaml'],fullfile(analysispath,'Bsn_link.yaml'));
                elseif contains(movepath,'Round2_right')
                    copyfile([movepath '.hdf5'],fullfile(analysispath,'CaV21_link.hdf5'));
                    copyfile([movepath '.yaml'],fullfile(analysispath,'CaV21_link.yaml'));
                end
            end % protein
            
            %%%%%%%% DBSCAN AND CLUSTERFILTER %%%%%%%%%%%%%%%%%%%%%%
            % grab files in analysis folder
            analysisfiles = dir(fullfile(analysispath,'*.hdf5'));
            analysispaths = fullfile({analysisfiles.folder}',{analysisfiles.name}');
            
            % DBSCAN and clusterFilter
            for db = 1:length(analysispaths)
                
                S = 'DBSCAN and clusterfilter on image %d of %d \n';
                fprintf(S,db,length(analysispaths))
                
                thistodbpath = analysispaths{db};
                if ~contains(thistodbpath,'CaV21')
                    
                    [dbpath,~,clusterpath] = picassoDbscanV057(thistodbpath, eps, minpts, pixelsize);
                    [in,out,dbs_hs] = clusterFilterV060(dbpath,clusterpath,...
                        std_range,'spotcheck',1,'minn',75);
                    
                    savefig(figure(1),[erase(dbpath,'.hdf5') '_CFgraph.fig'])
                    savefig(figure(2),[erase(dbpath,'.hdf5') '_CFspotcheck.fig'])
                    close all 
                    
                    RPTP(in,dbs_hs,erase(dbpath,'.hdf5'),'_CF_in');
                    RPTP(out,dbs_hs,erase(dbpath,'.hdf5'),'_CF_out');  
                    
                elseif contains(thistodbpath,'CaV21')
                    
                    [dbpath,~,clusterpath] = picassoDbscanV057(thistodbpath, eps, minpts, pixelsize);
                    [in,out,dbs_hs] = clusterFilterV060(dbpath,clusterpath,...
                        std_range,'spotcheck',1,'minn',20);
                    
                    savefig(figure(1),[erase(dbpath,'.hdf5') '_CFgraph.fig'])
                    savefig(figure(2),[erase(dbpath,'.hdf5') '_CFspotcheck.fig'])
                    close all
                    
                    RPTP(in,dbs_hs,erase(dbpath,'.hdf5'),'_CF_in');
                    RPTP(out,dbs_hs,erase(dbpath,'.hdf5'),'_CF_out');
                end
                
            end %dbscan/cf

            % make a final figure
            [~,CaV21,col] = Omniloader(fullfile(analysispath,'CaV21_link_dbscan_CF_in.hdf5'),'verbose',0);
            [~,Bsn,~] = Omniloader(fullfile(analysispath,'Bsn_link_dbscan_CF_in.hdf5'),'verbose',0);
            [~,PSD95,~] = Omniloader(fullfile(analysispath,'PSD95_link_dbscan_CF_in.hdf5'),'verbose',0);
            [~,Munc13,~] = Omniloader(fullfile(analysispath,'Munc13_link_dbscan_CF_in.hdf5'),'verbose',0);

            figure
            hold on
            scatter(PSD95(:,col.x),PSD95(:,col.y),1,'.b');
            scatter(Bsn(:,col.x),Bsn(:,col.y),1,'.g');
            scatter(Munc13(:,col.x),Munc13(:,col.y),1,'.c');
            scatter(CaV21(:,col.x),CaV21(:,col.y),1,'.m')
            axis equal
            legend({'PSD95','Bsn','Munc13','CaV21'})
            savefig(gcf,fullfile(analysispath,'allcolorsfinal.fig'))
            
            fprintf('Finished region %d \n',thisrgn)
            
        end %rgn
        
        fprintf('Finished imaging day %d \n',thisimday)
        
    end %imday
    
    fprintf('Finished week %d \n', thisweek)
    
end %week

fprintf('Processing is finished!')
