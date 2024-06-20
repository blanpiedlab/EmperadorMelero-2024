% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine

%% Load
load('FINAL_NC_Analysis_redoAdded.mat');
dataset = 'CaV';

%% Display fieldnames for synaptic data 
datatype = 'Syn';
subs = AllData.([dataset '_' datatype '_Dataset']);

names = fieldnames(subs);
v = [num2cell(1:size(names,1))' names];
disp(v)

%% Synapse areas (converted to nm^2)
datatype = 'Syn';
field = 'Channel_1_SynArea';
synarea_bsn = extractNCAnalyses(AllData, dataset, datatype, field).*160^2;

field = 'Channel_2_SynArea';
synarea_psd = extractNCAnalyses(AllData, dataset, datatype, field).*160^2;

field = 'Channel_3_SynArea';
synarea_munc = extractNCAnalyses(AllData, dataset, datatype, field).*160^2;

field = 'Channel_4_SynArea';
synarea_cav = extractNCAnalyses(AllData, dataset, datatype, field).*160^2;

aa_synareas = [synarea_munc, synarea_cav, synarea_bsn, synarea_psd];

%% Synapse autocorrelations
datatype = 'Syn';
field = 'Channel_1_AutoCorr';
ac_bsn = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_2_AutoCorr';
ac_psd = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_3_AutoCorr';
ac_munc = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_4_AutoCorr';
ac_cav = extractNCAnalyses(AllData, dataset, datatype, field);

aa_acs = [ac_munc(:,2:end); ac_cav(:,2:end); ac_bsn(:,2:end); ac_psd(:,2:end)]; % crop first value that we don't use

%% Synapse centrality sqrt (per loc)
datatype = 'Syn';
field = 'Channel_1_CentralitySqrt';
centsqrt_bsn = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_2_CentralitySqrt';
centsqrt_psd = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_3_CentralitySqrt';
centsqrt_munc = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_4_CentralitySqrt';
centsqrt_cav = extractNCAnalyses(AllData, dataset, datatype, field);

aa_centsqrts = nan(max([size(centsqrt_bsn,1),size(centsqrt_psd,1),size(centsqrt_munc,1),size(centsqrt_cav,1)]),4);
aa_centsqrts(1:size(centsqrt_munc,1),1) = centsqrt_munc;
aa_centsqrts(1:size(centsqrt_cav,1),2) = centsqrt_cav;
aa_centsqrts(1:size(centsqrt_bsn,1),3) = centsqrt_bsn;
aa_centsqrts(1:size(centsqrt_psd,1),4) = centsqrt_psd;

%% Synapse centrality sqrt (per synapse)
datatype = 'Syn';
field = 'Channel_1_LocNum';
locnum_bsn = extractNCAnalyses(AllData, dataset, datatype, field);
centsqrtpersyn_bsn = [];
start = 1;
for bsni = 1:size(locnum_bsn,1)
    stop = start + locnum_bsn(bsni) - 1;
    centsqrtpersyn_bsn = [centsqrtpersyn_bsn; mean(centsqrt_bsn(start:stop))];
    start = stop + 1;
end

field = 'Channel_2_LocNum';
locnum_psd = extractNCAnalyses(AllData, dataset, datatype, field);
centsqrtpersyn_psd = [];
start = 1;
for psdi = 1:size(locnum_psd,1)
    stop = start + locnum_psd(psdi) - 1;
    centsqrtpersyn_psd = [centsqrtpersyn_psd; mean(centsqrt_psd(start:stop))];
    start = stop + 1;
end

field = 'Channel_3_LocNum';
locnum_munc = extractNCAnalyses(AllData, dataset, datatype, field);
centsqrtpersyn_munc = [];
start = 1;
for munci = 1:size(locnum_munc,1)
    stop = start + locnum_munc(munci) - 1;
    centsqrtpersyn_munc = [centsqrtpersyn_munc; mean(centsqrt_munc(start:stop))];
    start = stop + 1;
end

field = 'Channel_4_LocNum';
locnum_cav = extractNCAnalyses(AllData, dataset, datatype, field);
centsqrtpersyn_cav = [];
start = 1;
for cavi = 1:size(locnum_cav,1)
    stop = start + locnum_cav(cavi) - 1;
    centsqrtpersyn_cav = [centsqrtpersyn_cav; mean(centsqrt_cav(start:stop))];
    start = stop + 1;
end

aa_centsqrtpersyns = [centsqrtpersyn_munc, centsqrtpersyn_cav, centsqrtpersyn_bsn, centsqrtpersyn_psd];
 
%% Crosscorrelations (cleaned; first column already removed)
datatype = 'Syn';
field = 'CrossCorr_Ch1_Ch2';
xc_bsnpsd = extractNCAnalyses(AllData, dataset, datatype, field);
xc_bsnpsd = cleanCrossCorr(xc_bsnpsd,1);

field = 'CrossCorr_Ch1_Ch3';
xc_bsnmunc = extractNCAnalyses(AllData, dataset, datatype, field);
xc_bsnmunc = cleanCrossCorr(xc_bsnmunc,1);

field = 'CrossCorr_Ch1_Ch4';
xc_bsncav = extractNCAnalyses(AllData, dataset, datatype, field);
xc_bsncav = cleanCrossCorr(xc_bsncav,1);

field = 'CrossCorr_Ch2_Ch3';
xc_psdmunc = extractNCAnalyses(AllData, dataset, datatype, field);
xc_psdmunc = cleanCrossCorr(xc_psdmunc,1);

field = 'CrossCorr_Ch2_Ch4';
xc_psdcav = extractNCAnalyses(AllData, dataset, datatype, field);
xc_psdcav = cleanCrossCorr(xc_psdcav,1);

field = 'CrossCorr_Ch3_Ch4';
xc_munccav = extractNCAnalyses(AllData, dataset, datatype, field);
xc_munccav = cleanCrossCorr(xc_munccav,1);

aa_xcs = [xc_bsnpsd; xc_bsnmunc; xc_bsncav; xc_psdmunc; xc_psdcav; xc_munccav];

%% Number of nanoclusters per synapse
datatype = 'Syn';
field = 'Channel_1_NCnum';
ncnum_bsn = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_2_NCnum';
ncnum_psd = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_3_NCnum';
ncnum_munc = extractNCAnalyses(AllData, dataset, datatype, field);

field = 'Channel_4_NCnum';
ncnum_cav = extractNCAnalyses(AllData, dataset, datatype, field);

aa_ncnum = [ncnum_munc, ncnum_cav, ncnum_bsn, ncnum_psd];

%% Display fieldnames for NC data
datatype = 'NC';
POI = 'Munc';

subs = AllData.([dataset '_' datatype '_Dataset']);
subsubs = subs.([POI '_Data']);

names = fieldnames(subsubs);
v = [num2cell(1:size(names,1))' names];
disp(v)

%% NC areas (convert to nm^2)
datatype = 'NC';
POI = 'Bsn';
field = [POI '_NC_Area'];
ncarea_bsn = extractNCAnalyses(AllData, dataset, datatype, field, POI).*160^2;

POI = 'PSD95';
field = [POI '_NC_Area'];
ncarea_psd = extractNCAnalyses(AllData, dataset, datatype, field, POI).*160^2;

POI = 'Munc';
field = [POI '_NC_Area'];
ncarea_munc = extractNCAnalyses(AllData, dataset, datatype, field, POI).*160^2;

POI = 'Cav';
field = [POI '_NC_Area'];
ncarea_cav = extractNCAnalyses(AllData, dataset, datatype, field, POI).*160^2;

aa_ncareas = nan(max([size(ncarea_bsn,1),size(ncarea_psd,1),size(ncarea_munc,1),size(ncarea_cav,1)]),4);
aa_ncareas(1:size(ncarea_munc,1),1) = ncarea_munc;
aa_ncareas(1:size(ncarea_cav,1),2) = ncarea_cav;
aa_ncareas(1:size(ncarea_bsn,1),3) = ncarea_bsn;
aa_ncareas(1:size(ncarea_psd,1),4) = ncarea_psd;

%% NC centrality
datatype = 'NC';
POI = 'Bsn';
field = [POI '_NC_EllipseCentralitySqrt_Real'];
nccentsqrt_bsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

POI = 'PSD95';
field = [POI '_NC_EllipseCentralitySqrt_Real'];
nccentsqrt_psd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

POI = 'Munc';
field = [POI '_NC_EllipseCentralitySqrt_Real'];
nccentsqrt_munc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

POI = 'Cav';
field = [POI '_NC_EllipseCentralitySqrt_Real'];
nccentsqrt_cav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

aa_nccentsqrts = nan(max([size(nccentsqrt_bsn,1),size(nccentsqrt_psd,1),size(nccentsqrt_munc,1),size(nccentsqrt_cav,1)]),4);
aa_nccentsqrts(1:size(nccentsqrt_munc,1),1) = nccentsqrt_munc;
aa_nccentsqrts(1:size(nccentsqrt_cav,1),2) = nccentsqrt_cav;
aa_nccentsqrts(1:size(nccentsqrt_bsn,1),3) = nccentsqrt_bsn;
aa_nccentsqrts(1:size(nccentsqrt_psd,1),4) = nccentsqrt_psd;

%% NC peak to peak
datatype = 'NC';
POI = 'Bsn';
field = [POI '_p2pnm_to_PSD95_Real'];
p2p_bsnpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Munc_Real'];
p2p_bsnmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Cav_Real'];
p2p_bsncav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

p2pbsn = [p2p_bsnpsd, p2p_bsnmunc, p2p_bsncav];


POI = 'PSD95';
field = [POI '_p2pnm_to_Bsn_Real'];
p2p_psdbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Munc_Real'];
p2p_psdmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Cav_Real'];
p2p_psdcav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

p2ppsd = [p2p_psdbsn, p2p_psdmunc, p2p_psdcav];


POI = 'Munc';
field = [POI '_p2pnm_to_Bsn_Real'];
p2p_muncbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_PSD95_Real'];
p2p_muncpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Cav_Real'];
p2p_munccav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

p2pmunc = [p2p_muncbsn, p2p_muncpsd, p2p_munccav];


POI = 'Cav';
field = [POI '_p2pnm_to_Bsn_Real'];
p2p_cavbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_PSD95_Real'];
p2p_cavpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_p2pnm_to_Munc_Real'];
p2p_cavmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

p2pcav = [p2p_cavbsn, p2p_cavpsd, p2p_cavmunc];


aa_p2ps = nan(max([size(p2pbsn,1),size(p2ppsd,1),size(p2pmunc,1),size(p2pcav,1)]),12);
aa_p2ps(1:size(p2pmunc,1),1:3) = p2pmunc;
aa_p2ps(1:size(p2pcav,1),4:6) = p2pcav;
aa_p2ps(1:size(p2pbsn,1),7:9) = p2pbsn;
aa_p2ps(1:size(p2ppsd,1),10:12) = p2ppsd;

%% Prep enrichments for outlier removal (first column cropped, trailing zero converted to nan)
datatype = 'NC';
POI = 'Bsn';
field = [POI '_Real_NC_Peak_Enrichment_to_Bsn_Dist'];
enrichprep_bsnbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_bsnbsn = prepEnrichmentForOutliers(enrichprep_bsnbsn,1);

field = [POI '_Real_NC_Peak_Enrichment_to_PSD95_Dist'];
enrichprep_bsnpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_bsnpsd = prepEnrichmentForOutliers(enrichprep_bsnpsd,1);

field = [POI '_Real_NC_Peak_Enrichment_to_Munc_Dist'];
enrichprep_bsnmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_bsnmunc = prepEnrichmentForOutliers(enrichprep_bsnmunc,1);

field = [POI '_Real_NC_Peak_Enrichment_to_Cav_Dist'];
enrichprep_bsncav = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_bsncav = prepEnrichmentForOutliers(enrichprep_bsncav,1);


POI = 'PSD95';
field = [POI '_Real_NC_Peak_Enrichment_to_Bsn_Dist'];
enrichprep_psdbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_psdbsn = prepEnrichmentForOutliers(enrichprep_psdbsn,1);

field = [POI '_Real_NC_Peak_Enrichment_to_PSD95_Dist'];
enrichprep_psdpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_psdpsd = prepEnrichmentForOutliers(enrichprep_psdpsd,1);

field = [POI '_Real_NC_Peak_Enrichment_to_Munc_Dist'];
enrichprep_psdmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_psdmunc = prepEnrichmentForOutliers(enrichprep_psdmunc,1);

field = [POI '_Real_NC_Peak_Enrichment_to_Cav_Dist'];
enrichprep_psdcav = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_psdcav = prepEnrichmentForOutliers(enrichprep_psdcav,1);


POI = 'Munc';
field = [POI '_Real_Peak_NC_Enrichment_to_Bsn_Dist']; % order of nc/peak changed but its the same field
enrichprep_muncbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_muncbsn = prepEnrichmentForOutliers(enrichprep_muncbsn,1);

field = [POI '_Real_Peak_NC_Enrichment_to_PSD95_Dist'];
enrichprep_muncpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_muncpsd = prepEnrichmentForOutliers(enrichprep_muncpsd,1);

field = [POI '_Real_Peak_NC_Enrichment_to_Munc_Dist'];
enrichprep_muncmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_muncmunc = prepEnrichmentForOutliers(enrichprep_muncmunc,1);

field = [POI '_Real_Peak_NC_Enrichment_to_Cav_Dist'];
enrichprep_munccav = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_munccav = prepEnrichmentForOutliers(enrichprep_munccav,1);


POI = 'Cav';
field = [POI '_Real_Peak_NC_Enrichment_to_Bsn_Dist'];
enrichprep_cavbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_cavbsn = prepEnrichmentForOutliers(enrichprep_cavbsn,1);

field = [POI '_Real_Peak_NC_Enrichment_to_PSD95_Dist'];
enrichprep_cavpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_cavpsd = prepEnrichmentForOutliers(enrichprep_cavpsd,1);

field = [POI '_Real_Peak_NC_Enrichment_to_Munc_Dist'];
enrichprep_cavmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_cavmunc = prepEnrichmentForOutliers(enrichprep_cavmunc,1);

field = [POI '_Real_Peak_NC_Enrichment_to_Cav_Dist'];
enrichprep_cavcav = extractNCAnalyses(AllData, dataset, datatype, field, POI);
enrichprep_cavcav = prepEnrichmentForOutliers(enrichprep_cavcav,1);

%% Final enrichments and EIs (outliers determined in PRISM with ROUT, pasted into cutoffs as row vector)
cutoffs = [3.67948	2.89855	2.01737	1.96945	2.10191	2.27704	2.13483	2.25346	2.10472	2.24678	2.41033	2.16047	2.12511	2.20538	2.18778	2.15311	2.22222	2.14286	2.46773	2.44399	2.57143	2.38095	2.37154	2.37718	2.70270	2.64609	2.67857	2.30179	3.03030	3.33333	2.90557	3.30806];
enrich_bsnbsn = smoothEnrichment(enrichprep_bsnbsn,cutoffs);
%ei_bsnbsn = mean(enrich_bsnbsn(:,1:5),2,'omitnan');

cutoffs = [3.00000	3.27238	3.04569	3.01109	2.79232	2.78364	2.72290	2.51917	2.74949	2.83353	2.57701	2.92758	3.02714	2.95463	2.99667	2.92050	3.03030	3.18021	3.08008	3.08300	3.36574	3.73250	3.81356	3.82746	3.89610	3.54275	4.34783	4.53686	5.55556	4.13223	5.34351	8.45070];
enrich_bsnpsd = smoothEnrichment(enrichprep_bsnpsd,cutoffs);
ei_bsnpsd = mean(enrich_bsnpsd(:,1:5),2,'omitnan');

cutoffs = [4.40868	4.35631	3.95991	3.68852	3.15851	3.11938	2.89140	2.92065	3.12306	2.71084	2.63158	2.64094	3.08422	2.84765	2.67631	2.82759	2.71669	2.90387	3.01766	2.42329	3.00752	2.51828	2.55134	2.42915	2.61122	2.49792	2.23325	2.44565	2.59965	2.13282	2.43463	3.03084];
enrich_bsnmunc = smoothEnrichment(enrichprep_bsnmunc,cutoffs);
ei_bsnmunc = mean(enrich_bsnmunc(:,1:5),2,'omitnan');

cutoffs = [3.46535	4.32692	4.73422	4.25197	4.32990	3.77020	3.70968	3.88170	3.98734	3.83292	3.49345	3.71901	3.79747	3.86740	3.11688	3.50195	3.39806	3.30579	3.42857	3.18302	3.30954	3.19829	3.27869	2.83912	3.34928	3.54610	3.06122	3.20856	2.92805	3.34572	3.73585	2.42616];
enrich_bsncav = smoothEnrichment(enrichprep_bsncav,cutoffs);
ei_bsncav = mean(enrich_bsncav(:,1:5),2,'omitnan');

aa_enrichbsns = [enrich_bsnbsn; enrich_bsnpsd; enrich_bsnmunc; enrich_bsncav];
aa_eis_bsn = [ei_bsnpsd, ei_bsnmunc, ei_bsncav];


cutoffs = [4.35931	3.58337	2.94915	2.91126	2.59827	2.66459	2.72134	2.53989	2.72057	2.70068	2.57813	2.88462	2.84464	2.85132	2.51046	2.86598	2.99472	2.48447	3.13808	2.64145	3.34780	3.73415	2.90404	3.46498	3.44828	3.19761	3.48837	2.48246	3.58065	3.04569	2.50256	3.60577];
enrich_psdbsn = smoothEnrichment(enrichprep_psdbsn,cutoffs);
ei_psdbsn = mean(enrich_psdbsn(:,1:5),2,'omitnan');

cutoffs = [4.56885	3.04816	2.35195	2.27743	2.13946	2.34375	2.20723	2.07357	2.38474	2.38095	2.24719	2.10186	2.18571	2.16606	2.01027	2.14550	2.08333	2.41287	2.07565	2.29358	2.12766	2.26415	2.15827	2.27273	2.52101	2.29008	2.09581	2.15827	2.57143	2.72480	2.72727	2.02703];
enrich_psdpsd = smoothEnrichment(enrichprep_psdpsd,cutoffs);
%ei_psdpsd = mean(enrich_psdpsd(:,1:5),2,'omitnan');

cutoffs = [4.69028	5.07302	4.01338	4.13793	3.68224	3.34169	3.20525	3.41081	2.84939	3.14361	3.21470	2.94530	2.90537	2.88462	3.23213	2.83067	2.99003	3.56779	3.06905	2.61745	2.93958	3.19878	2.88462	2.46385	3.25123	2.50232	2.40320	3.29670	1.78268	2.67431	2.29520	1.65837];
enrich_psdmunc = smoothEnrichment(enrichprep_psdmunc,cutoffs);
ei_psdmunc = mean(enrich_psdmunc(:,1:5),2,'omitnan');

cutoffs = [3.66748	4.55531	4.54122	4.27566	4.10959	3.34728	3.52273	3.38674	3.70722	3.41834	2.87411	3.63541	3.65854	3.44828	3.80952	3.65854	3.26087	3.48162	3.75671	3.92523	3.38983	3.25301	3.00214	3.44828	3.54331	3.70066	3.57798	3.29068	3.08087	3.64217	2.84834	2.59570];
enrich_psdcav = smoothEnrichment(enrichprep_psdcav,cutoffs);
ei_psdcav = mean(enrich_psdcav(:,1:5),2,'omitnan');

aa_enrichpsds = [enrich_psdbsn; enrich_psdpsd; enrich_psdmunc; enrich_psdcav];
aa_eis_psd = [ei_psdbsn, ei_psdmunc, ei_psdcav];


cutoffs = [3.52505	2.86156	2.29913	2.26907	2.22140	2.08556	2.41935	1.98226	2.51908	2.16082	2.38754	2.32558	2.08801	2.28311	2.19468	2.25466	2.45902	2.36220	2.71903	2.58621	2.51497	2.72727	2.54669	2.59516	2.33766	2.22987	2.67857	2.94118	2.83019	2.41854	3.10881	2.34375];
enrich_muncbsn = smoothEnrichment(enrichprep_muncbsn,cutoffs);
ei_muncbsn = mean(enrich_muncbsn(:,1:5),2,'omitnan');

cutoffs = [4.36508	3.38770	3.35196	3.22027	2.85326	3.00222	2.62887	2.69352	2.56035	2.44051	2.61673	2.25564	2.51601	2.86244	2.72727	3.03030	3.00200	3.26087	3.03468	4.05405	3.93586	4.22886	4.46927	4.19580	4.78088	4.41989	4.17311	4.26136	7.31707	7.32601	5.09491	4.51128];
enrich_muncpsd = smoothEnrichment(enrichprep_muncpsd,cutoffs);
ei_muncpsd = mean(enrich_muncpsd(:,1:5),2,'omitnan');

cutoffs = [4.84282	3.93701	3.51491	3.10581	3.05060	2.69416	2.52746	2.40783	2.58993	2.31327	2.23287	2.31130	2.38241	2.33909	2.26827	2.27812	2.66667	2.01632	2.72577	2.52996	2.12371	1.96137	1.80389	2.50569	1.78674	2.41228	2.32420	2.40000	1.95713	1.60234	2.01806	1.91975];
enrich_muncmunc = smoothEnrichment(enrichprep_muncmunc,cutoffs);
%ei_muncmunc = mean(enrich_muncmunc(:,1:5),2,'omitnan');

cutoffs = [5.33333	5.17928	4.50867	4.07932	3.82746	3.41207	3.64583	3.43035	3.46420	3.36323	3.00375	3.25581	3.71367	3.37942	3.22896	3.01003	3.16056	3.56612	3.02521	3.03933	3.43277	3.52941	3.34572	3.11111	3.33333	2.61166	3.26087	3.02077	3.22581	2.80851	2.91262	2.56959];
enrich_munccav = smoothEnrichment(enrichprep_munccav,cutoffs);
ei_munccav = mean(enrich_munccav(:,1:5),2,'omitnan');

aa_enrichmuncs = [enrich_muncbsn; enrich_muncpsd; enrich_muncmunc; enrich_munccav];
aa_eis_munc = [ei_muncbsn, ei_muncpsd, ei_munccav];


cutoffs = [3.39623	3.11419	2.69830	2.24138	2.36981	2.44565	2.13704	1.92373	2.23547	2.52101	1.95298	2.04280	2.11885	2.28517	2.56410	2.60870	2.03666	2.60870	2.40642	2.94118	2.98211	2.32558	2.70270	2.26415	2.64901	2.08968	2.72727	2.79720	2.61324	2.10379	3.48837	2.24426];
enrich_cavbsn = smoothEnrichment(enrichprep_cavbsn,cutoffs);
ei_cavbsn = mean(enrich_cavbsn(:,1:5),2,'omitnan');

cutoffs = [2.93734	3.35526	2.96296	2.88000	2.77980	2.62697	2.48157	2.59638	2.41503	2.25847	2.73764	2.38411	2.75574	2.77778	2.71710	2.54561	2.89575	3.35508	3.45061	3.39175	3.21888	3.52719	2.41935	3.52941	3.50028	4.87692	6.48199	7.75717	6.00000	6.25000	4.79651	6.12245];
enrich_cavpsd = smoothEnrichment(enrichprep_cavpsd,cutoffs);
ei_cavpsd = mean(enrich_cavpsd(:,1:5),2,'omitnan');

cutoffs = [4.59184	3.69863	3.37079	2.75385	3.01653	2.60816	2.73816	2.69836	2.13662	1.98939	2.16426	2.45399	2.42188	2.45148	2.04380	3.14075	2.85395	2.76400	2.58621	2.52606	2.17391	2.52487	2.63158	1.79356	2.12332	2.17658	2.12353	1.98407	2.19597	2.68520	2.38938	2.30216];
enrich_cavmunc = smoothEnrichment(enrichprep_cavmunc,cutoffs);
ei_cavmunc = mean(enrich_cavmunc(:,1:5),2,'omitnan');

cutoffs = [10.0159	5.10949	3.54610	3.43075	2.89700	3.21970	2.47934	2.47841	2.26415	2.66075	2.19713	2.65487	2.34681	2.56024	2.48242	2.76833	2.79070	2.59034	2.90248	2.50397	2.93638	2.83480	2.82220	3.32226	1.91599	2.12202	2.65487	2.34719	2.37288	2.55910	2.32996	1.73675];
enrich_cavcav = smoothEnrichment(enrichprep_cavcav,cutoffs);
%ei_cavcav = mean(enrich_cavcav(:,1:5),2,'omitnan');

aa_enrichcavs = [enrich_cavbsn; enrich_cavpsd; enrich_cavmunc; enrich_cavcav];
aa_eis_cav = [ei_cavbsn, ei_cavpsd, ei_cavmunc];

aa_eis = nan(max([size(aa_eis_bsn,1),size(aa_eis_psd,1),size(aa_eis_munc,1),size(aa_eis_cav,1)]),12);
aa_eis(1:size(aa_eis_munc,1),1:3) = aa_eis_munc;
aa_eis(1:size(aa_eis_cav,1),4:6) = aa_eis_cav;
aa_eis(1:size(aa_eis_bsn,1),7:9) = aa_eis_bsn;
aa_eis(1:size(aa_eis_psd,1),10:12) = aa_eis_psd;

%% mean isAligned measures
datatype = 'NC';
POI = 'Bsn';
field = [POI '_NC_Peak_Mean_isAligned_to_PSD95_Dist'];
isaligned_bsnpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Munc_Dist'];
isaligned_bsnmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Cav_Dist'];
isaligned_bsncav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

isalignedbsn = [isaligned_bsnpsd, isaligned_bsnmunc, isaligned_bsncav];


POI = 'PSD95';
field = [POI '_NC_Peak_Mean_isAligned_to_Bsn_Dist'];
isaligned_psdbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Munc_Dist'];
isaligned_psdmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Cav_Dist'];
isaligned_psdcav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

isalignedpsd = [isaligned_psdbsn, isaligned_psdmunc, isaligned_psdcav];


POI = 'Munc';
field = [POI '_NC_Peak_Mean_isAligned_to_Bsn_Dist'];
isaligned_muncbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_PSD95_Dist'];
isaligned_muncpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Cav_Dist'];
isaligned_munccav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

isalignedmunc = [isaligned_muncbsn, isaligned_muncpsd, isaligned_munccav];


POI = 'Cav';
field = [POI '_NC_Peak_Mean_isAligned_to_Bsn_Dist'];
isaligned_cavbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_PSD95_Dist'];
isaligned_cavpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_NC_Peak_Mean_isAligned_to_Munc_Dist'];
isaligned_cavmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

isalignedcav = [isaligned_cavbsn, isaligned_cavpsd, isaligned_cavmunc];


aa_isaligneds = nan(max([size(isalignedbsn,1),size(isalignedpsd,1),size(isalignedmunc,1),size(isalignedcav,1)]),12);
aa_isaligneds(1:size(isalignedmunc,1),1:3) = isalignedmunc;
aa_isaligneds(1:size(isalignedcav,1),4:6) = isalignedcav;
aa_isaligneds(1:size(isalignedbsn,1),7:9) = isalignedbsn;
aa_isaligneds(1:size(isalignedpsd,1),10:12) = isalignedpsd;

%% AI
datatype = 'NC';
POI = 'Bsn';
field = [POI '_AI_to_PSD95_Real'];
ai_bsnpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Munc_Real'];
ai_bsnmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Cav_Real'];
ai_bsncav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

aibsn = [ai_bsnpsd, ai_bsnmunc, ai_bsncav];


POI = 'PSD95';
field = [POI '_AI_to_Bsn_Real'];
ai_psdbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Munc_Real'];
ai_psdmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Cav_Real'];
ai_psdcav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

aipsd = [ai_psdbsn, ai_psdmunc, ai_psdcav];


POI = 'Munc';
field = [POI '_AI_to_Bsn_Real'];
ai_muncbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_PSD95_Real'];
ai_muncpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Cav_Real'];
ai_munccav = extractNCAnalyses(AllData, dataset, datatype, field, POI);

aimunc = [ai_muncbsn, ai_muncpsd, ai_munccav];


POI = 'Cav';
field = [POI '_AI_to_Bsn_Real'];
ai_cavbsn = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_PSD95_Real'];
ai_cavpsd = extractNCAnalyses(AllData, dataset, datatype, field, POI);

field = [POI '_AI_to_Munc_Real'];
ai_cavmunc = extractNCAnalyses(AllData, dataset, datatype, field, POI);

aicav = [ai_cavbsn, ai_cavpsd, ai_cavmunc];


aa_ai = nan(max([size(aibsn,1),size(aipsd,1),size(aimunc,1),size(aicav,1)]),12);
aa_ai(1:size(aimunc,1),1:3) = aimunc;
aa_ai(1:size(aicav,1),4:6) = aicav;
aa_ai(1:size(aibsn,1),7:9) = aibsn;
aa_ai(1:size(aipsd,1),10:12) = aipsd;

%% Some specific measurements - isaligned %s, enrichments coded by various isaligned
% enrichments coded by isaligned
enrich_munccav_by_muncpsd_aligned = enrich_munccav;
enrich_munccav_by_muncpsd_notaligned = enrich_munccav;
enrich_munccav_by_muncpsd_undetermined = enrich_munccav;
enrich_munccav_by_muncpsd_aligned(isaligned_muncpsd~=1,:) = nan;
enrich_munccav_by_muncpsd_notaligned(isaligned_muncpsd~=-1,:) = nan;
enrich_munccav_by_muncpsd_undetermined(isaligned_muncpsd~=0,:) = nan;
aa_enrich_munccav_by_muncpsd_isAligned = [enrich_munccav_by_muncpsd_aligned; enrich_munccav_by_muncpsd_notaligned; enrich_munccav_by_muncpsd_undetermined];

enrich_munccav_by_munccav_aligned = enrich_munccav;
enrich_munccav_by_munccav_notaligned = enrich_munccav;
enrich_munccav_by_munccav_undetermined = enrich_munccav;
enrich_munccav_by_munccav_aligned(isaligned_munccav~=1,:) = nan;
enrich_munccav_by_munccav_notaligned(isaligned_munccav~=-1,:) = nan;
enrich_munccav_by_munccav_undetermined(isaligned_munccav~=0,:) = nan;
aa_enrich_munccav_by_munccav_isAligned = [enrich_munccav_by_munccav_aligned; enrich_munccav_by_munccav_notaligned; enrich_munccav_by_munccav_undetermined];

enrich_muncpsd_by_munccav_aligned = enrich_muncpsd;
enrich_muncpsd_by_munccav_notaligned = enrich_muncpsd;
enrich_muncpsd_by_munccav_undetermined = enrich_muncpsd;
enrich_muncpsd_by_munccav_aligned(isaligned_munccav~=1,:) = nan;
enrich_muncpsd_by_munccav_notaligned(isaligned_munccav~=-1,:) = nan;
enrich_muncpsd_by_munccav_undetermined(isaligned_munccav~=0,:) = nan;
aa_enrich_muncpsd_by_munccav_isAligned = [enrich_muncpsd_by_munccav_aligned; enrich_muncpsd_by_munccav_notaligned; enrich_muncpsd_by_munccav_undetermined];

load('enrich_muncpsd.mat');
load('aa_isaligneds.mat');
enrich_muncpsd_by_muncpsd_aligned = enrich_muncpsd;
enrich_muncpsd_by_muncpsd_notaligned = enrich_muncpsd;
enrich_muncpsd_by_muncpsd_undetermined = enrich_muncpsd;
enrich_muncpsd_by_muncpsd_aligned(aa_isaligneds(1:350,2)~=1,:) = nan;
enrich_muncpsd_by_muncpsd_notaligned(aa_isaligneds(1:350,2)~=-1,:) = nan;
enrich_muncpsd_by_muncpsd_undetermined(aa_isaligneds(1:350,2)~=0,:) = nan;
aa_enrich_muncpsd_by_muncpsd_isAligned = [enrich_muncpsd_by_muncpsd_aligned; enrich_muncpsd_by_muncpsd_notaligned; enrich_muncpsd_by_muncpsd_undetermined];

% isaligned percents of total
n = size(isaligned_munccav,1);
aa_isaligned_munccav_pcts = [sum(isaligned_munccav==1)/n; sum(isaligned_munccav==-1)/n; sum(isaligned_munccav==0)/n];

n = size(aa_isaligneds(1:350,2),1);
aa_isaligned_muncpsd_pcts = [sum(aa_isaligneds(1:350,2)==1)/n; sum(aa_isaligneds(1:350,2)==-1)/n; sum(aa_isaligneds(1:350,2)==0)/n];


% scatter plots
aa_ei_munccav_by_ei_muncpsd = [ei_munccav, ei_muncpsd];
aa_p2p_munccav_by_ei_muncpsd = [p2p_munccav, ei_muncpsd];
aa_ei_muncpsd_by_p2p_munccav = [ei_muncpsd, p2p_munccav];
aa_p2p_muncpsd_by_ei_munccav = [p2p_muncpsd, ei_munccav];
aa_p2p_muncpsd_by_p2p_munccav = [p2p_muncpsd, p2p_munccav];

%% Save data
varsinwksp = who;
for i = 1:numel(varsinwksp)
    
    thisvar = varsinwksp{i};
    if startsWith(thisvar,'aa_') || startsWith(thisvar,'enrich_') || startsWith(thisvar,'enrichprep_')
        save([thisvar '.mat'],thisvar)
    end
    
end
