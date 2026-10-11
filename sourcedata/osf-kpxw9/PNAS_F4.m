function PNAS_F4(datadir)
% PNAS_F4(datadir)
% code to produce Figure F4
%
% datadir ... path to 'sharing' folder (e.g., 'C:\mydata\sharing')
%
% note: script is based on fieldtrip-20170206
%

close all

if contains(path,'fieldtrip')
    [~, ft_dir] = ft_version;
    ft_defaults;
else
    error('FieldTrip not in path');
end

addpath(fullfile(datadir,'auxiliary'));
load(fullfile(datadir,'auxiliary','mycoldmap'),'mycoldmap')

%% settings

conditionnames      = {'Remember' 'Forgot'};
subjects            = {'s01' 's02' 's03' 's04' 's05' 's06' 's07' 's08' 's09' 's10' 's11' 's12' 's13' 's14' 's15' 's16' 's17' 's18' 's19' 's20'};
TOI                 = [-.5 1.6]; % global interval of interest
ubermat             = [];

%% enter subjects loop

for isubject = 1:numel(subjects)
        
    fprintf('working on subject %02d/%02d\n',isubject,numel(subjects));

    % Recall data
    tmp     = load(fullfile(datadir,'EEG',sprintf('s%02d_recall',isubject)));
    data    = tmp.data;
    clear tmp
    
    % reduce data size
    cfg             = [];
    cfg.precision   = 'single';
    data            = ft_preprocessing(cfg, data);

    % reduce data to time of interest plus 1 s buffer
    cfg             = [];
    cfg.latency     = [TOI(1)-1, TOI(2)+1];
    data            = ft_selectdata(cfg, data);
    
    % bring into matrix format
    cfg             = [];
    cfg.removemean  = 'no';
    cfg.keeptrials  = 'yes';
    data            = ft_timelockanalysis(cfg, data);
    
    %% go to source space
    
    fprintf('Moving to source space\n')
    
    % Load beamformer weights
    tmp      = load(fullfile(datadir,'EEG',sprintf('s%02d_beamformer_FEM',isubject)));
    W        = tmp.W;
    nweights = size(W,1);
    clear tmp

    trial = zeros(size(data.trial,1), nweights, size(data.trial,3));
    
    for ii=1:size(data.trial,1)
        trial(ii,:,:) = W * squeeze(data.trial(ii,:,:));
    end
    data.trial   = trial;
    clear trial W
    
    data.label = arrayfun( @(x) ['s' num2str(x)], 1:nweights, 'Un', 0);
    
    if isfield(data,'elec')
        data = rmfield(data,'elec');
    end
    
    %% frequency analysis using Hanning tapers
    
    fsample = 1/mode(diff(data.time));
    
    % Indices of desired time interval
    toi_idx  = nearest(data.time,TOI(1)):round(.050*fsample):nearest(data.time,TOI(2));
    
    cfg              = [];
    cfg.channel      = 'all';
    cfg.output       = 'pow';
    cfg.pad          = 'nextpow2';
    cfg.keeptrials   = 'yes';   
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 8:12;
    cfg.t_ftimwin    = 5./cfg.foi;
    cfg.toi          = data.time(toi_idx);
    freq            = ft_freqanalysis(cfg, data);
    clear data
    
    %% trimmean and z score
        
    trlinfo  = cell2mat(freq.trialinfo);
    Response = [trlinfo.Response];
    clear trlinfo

    for icondition = 1:numel(conditionnames)
        switch icondition
            case 1 % 'Remember'
                trialselection = ismember(Response,'left');
            case 2 % 'Forgot'
                trialselection = ~ismember(Response,'left');
        end
        
        dat = freq.powspctrm(trialselection==1,:,:,:);
        
        thresh          = prctile(dat,90,1);
        thresh          = repmat(thresh,[sum(trialselection) 1 1 1]);
        outl            = dat > thresh;
        clear thresh
        dat(outl==1)    = nan;
        clear outl            
        
        freq.powspctrm(trialselection==1,:,:,:) = dat;
        
        clear dat
    end
    
    %--- z score 
    pow  = freq.powspctrm;
    powz = nan(size(pow));
    
    for ichan = 1:size(pow,2)
        for ifreq = 1:size(pow,3)
            d = squeeze(pow(:,ichan,ifreq,:));
            m = nanmean(d(:));
            s = nanstd(d(:));
            powz(:,ichan,ifreq,:) = (pow(:,ichan,ifreq,:) - m)./s;
            clear d m s
        end
    end
    clear pow
    freq.powspctrm = powz;
    clear powz
    
    %% condition loop
    
    for icondition = 1:numel(conditionnames)     
        switch icondition
            case 1 % 'Remember'
                trialselection = ismember(Response,'left');
            case 2 % 'Forgot'
                trialselection = ~ismember(Response,'left');
        end
        
        cfg         = [];
        cfg.trials  = find(trialselection);
        cond        = ft_selectdata(cfg, freq);
        
        cond = rmfield(cond, 'trialinfo');
        
        cfg                 = [];
        cfg.nanmean         = 'yes';
        cfg.avgoverfreq     = 'yes';
        cfg.avgoverrpt      = 'yes';
        cfg.avgovertime     = 'no';
        cond                = ft_selectdata(cfg, cond);
        
        ubermat{isubject}{icondition} = cond;
        
        clear cond
    end
    clear freq
end

timeaxis = ubermat{1}{1}.time;

%% prepare source analyses and plotting 

load(fullfile(datadir,'EEG','template_source'),'template_source');

template_mri        = ft_read_mri([ft_dir '/template/headmodel/standard_mri.mat']);
templatedir         = fullfile(ft_dir, 'template', 'sourcemodel');
template_grid       = load(fullfile(templatedir, 'standard_sourcemodel3d10mm'));
template_grid       = template_grid.sourcemodel;
template_source.dim = template_grid.dim;
template_source.pos = template_grid.pos;

%% Figure 4A: grand average from 0-1.5 s

source_TOI  = [0 1.5];
surfaces    = {'surface_pial_left.mat'     'surface_pial_left.mat'     'surface_pial_both.mat'};
inflates    = {'surface_inflated_left.mat' 'surface_inflated_left.mat' 'surface_inflated_both.mat'};

sources1 = cell(numel(subjects),1);
sources2 = cell(numel(subjects),1);

for isubject = 1:numel(subjects)
    
    tmp                   = template_source;
    tmp.time              = timeaxis;
    tmp.pow               = zeros(numel(tmp.inside),numel(timeaxis));
    tmp.pow(tmp.inside,:) = squeeze(ubermat{isubject}{1}.powspctrm);
    sources1{isubject,1}  = tmp;
    tmp.pow               = zeros(numel(tmp.inside),numel(timeaxis));
    tmp.pow(tmp.inside,:) = squeeze(ubermat{isubject}{2}.powspctrm);
    sources2{isubject,1}  = tmp;
    clear tmp
end

cfg                     = [];
cfg.latency             = source_TOI;
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.statistic           = 'depsamplesT';
cfg.parameter           = 'pow';
cfg.numrandomization    = 1000;
cfg.alpha               = 0.05;
cfg.correcttail         = 'alpha';
cfg.clusteralpha        = 0.05;
cfg.tail                = 0;

cfg.design                         = [];
cfg.design(1,1:2*numel(subjects))  = [ones(1,numel(subjects)) 2*ones(1,numel(subjects))];
cfg.design(2,1:2*numel(subjects))  = [1:numel(subjects) 1:numel(subjects)];

cfg.ivar                = 1;
cfg.uvar                = 2;

sourcestat = ft_sourcestatistics(cfg, sources1{:,1}, sources2{:,1});

mask = zeros(size(sourcestat.mask));
mask(sourcestat.mask==1) = 1;

mystat          = sourcestat;
mystat.stat     = sum((sourcestat.stat .* (sourcestat.negclusterslabelmat ==1)) .* mask,2);

cfg             = [];
cfg.parameter   = 'stat';
int_source      = ft_sourceinterpolate(cfg, mystat, template_mri);

for iview = 1:numel(surfaces)
    
    cfg                 = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'stat';
    cfg.funcolorlim     = 'minzero';
    cfg.maskparameter   = cfg.funparameter;
    cfg.maskstyle       = 'opacity';
    cfg.opacitymap      = 'rampdown';
    cfg.projmethod      = 'project';
    cfg.surffile        = surfaces{iview};
    cfg.surfinflated    = inflates{iview};
    cfg.camlight        = 'no';
    cfg.funcolormap     = flipud(mycoldmap);
    
    ft_sourceplot(cfg, int_source);
    
    switch iview
        case 1
            view(90,0)
        case 2
            view(-90,0)
        case 3
            view(0,-90)
    end
    h = colorbar;
    ylabel(h, 't value')
    axis off
    axis equal;
    set(gca,'fontsize',8)
    caxis([-30 0])
end

%% Figure 4C: multiple ROI effect size time series

ROInames = {'HIPP','PPC','vmPFC'};

ROI        = {};
ROI{1}     = {...
    'Hippocampus_L'  'Hippocampus_R'};
ROI{2}     = {...
    'Precuneus_L' 'Cingulum_Post_L' 'Precuneus_R' 'Cingulum_Post_R'};
ROI{3}     = {...
    'Rectus_L'  'Frontal_Sup_Orb_L' 'Frontal_Med_Orb_L' 'Rectus_R'  'Frontal_Sup_Orb_R' 'Frontal_Med_Orb_R'}; %

ROI_TOI = [-.5 1.6];
atlas   = ft_read_atlas([ft_dir '/template/atlas/aal/ROI_MNI_V4.nii']);

timeaxis    = ubermat{1}{1}.time;
toi_idx     = nearest(timeaxis, ROI_TOI(1)):nearest(timeaxis, ROI_TOI(2));
timeaxis    = timeaxis(toi_idx);
groupdat    = {};

for isubject=1:numel(subjects)
    
    for icondition=1:numel(conditionnames)
        
        placeholder = zeros(numel(template_source.inside),numel(toi_idx));
        placeholder(template_source.inside,:) = squeeze(ubermat{isubject}{icondition}.powspctrm(:,:,toi_idx));
        groupdat{isubject,icondition}  = placeholder;        
    end
end

%--- get voxels for ROI
roi_grid   = ft_convert_units(template_grid, 'mm');

c = {};
for iroi = 1:numel(ROI)
    % find the ROI
    cfg             = [];
    cfg.atlas       = atlas;
    cfg.inputcoord  = 'mni';
    cfg.roi         = ROI{iroi};
    mask            = ft_volumelookup(cfg, roi_grid);
    
    % get data
    roi_idx = mask(:);
    roi_idx = roi_idx .* any(sourcestat.mask,2);

    fprintf('%d voxels in ROI\n',sum(roi_idx))
    
    for isubject = 1:size(groupdat,1)
        c{isubject,iroi,1}          = struct;
        c{isubject,iroi,1}.dimord   = 'rpt_chan_time';
        c{isubject,iroi,1}.label    = {'ROI'};
        c{isubject,iroi,1}.time     = timeaxis;
        c{isubject,iroi,1}.trial    = mean(groupdat{isubject,1}(roi_idx==1,:));
    end
    
    for isubject = 1:size(groupdat,1)
        c{isubject,iroi,2}          = struct;
        c{isubject,iroi,2}.dimord   = 'rpt_chan_time';
        c{isubject,iroi,2}.label    = {'ROI'};
        c{isubject,iroi,2}.time     = timeaxis;
        c{isubject,iroi,2}.trial    = mean(groupdat{isubject,2}(roi_idx==1,:));
    end
end

%--- plot
dat = [];
for iroi = 1:numel(ROI)
    for isubject=1:numel(subjects)
        dat(isubject,:,iroi) = c{isubject,iroi,1}.trial-c{isubject,iroi,2}.trial;
    end
end

[~,~,~,stat] = ttest(dat);

tvals = squeeze(stat.tstat);
tcrit = tinv(.05,mode(stat.df(:)));

regioncolors = {'g','b','r'};
figure;
hold on
l = [];
for ireg=1:size(tvals,2)
    l(ireg)=plot(timeaxis,tvals(:,ireg),'color',regioncolors{ireg},'linewidth',3);
end

hlinedashed(tcrit)
vline(0)
hline(0)
xlim([.5 1.5])
legend(l,ROInames,'location','best')
set(gca,'FontSize',8,'FontName','Helvetica','TickDir','out','Box','off')
xlabel('time (s)')
ylabel('t vals')

%% Figure 4B: multiple time windows from 100 ms - 1600 ms in 100 ms steps

source_TOIs = [0:.100:1.5;.100:.100:1.6]';

for iTOI = 1:size(source_TOIs,1)
    
    source_TOI_idx  = nearest(timeaxis, source_TOIs(iTOI,1)):nearest(timeaxis, source_TOIs(iTOI,2));
        
    sources1 = cell(numel(subjects),1);
    sources2 = cell(numel(subjects),1);
    
    for isubject = 1:numel(subjects)
        
        tmp                 = template_source;
        tmp.pow             = zeros(numel(tmp.inside),1);
        tmp.pow(tmp.inside) = squeeze(mean(ubermat{isubject}{1}.powspctrm(:,:,source_TOI_idx),3));
        
        sources1{isubject,1} = tmp;
        
        tmp.pow             = zeros(numel(tmp.inside),1);
        tmp.pow(tmp.inside) = squeeze(mean(ubermat{isubject}{2}.powspctrm(:,:,source_TOI_idx),3));
        
        sources2{isubject,1} = tmp;
        
        clear tmp
    end
    
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.statistic           = 'depsamplesT';
    cfg.parameter           = 'pow';
    cfg.numrandomization    = 1000;
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.clusteralpha        = 0.05;
    cfg.tail                = 0;
    cfg.ivar                = 1; 
    cfg.uvar                = 2;
    
    cfg.design                         = [];
    cfg.design(1,1:2*numel(subjects))  = [ones(1,numel(subjects)) 2*ones(1,numel(subjects))];
    cfg.design(2,1:2*numel(subjects))  = [1:numel(subjects) 1:numel(subjects)];
    
    sourcestat = ft_sourcestatistics(cfg, sources1{:,1}, sources2{:,1});
        
    if any(sourcestat.mask)

    mystat              = sourcestat;
    mystat.stat         = sourcestat.stat .* (sourcestat.negclusterslabelmat == 1);
    
    cfg                 = [];
    cfg.parameter       = 'stat';
    int_source          = ft_sourceinterpolate(cfg, mystat, template_mri);
    
    cfg                 = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'stat';
    cfg.funcolorlim     = 'minzero';
    cfg.maskparameter   = cfg.funparameter;
    cfg.maskstyle       = 'opacity';
    cfg.opacitymap      = 'rampdown';
    cfg.projmethod      = 'project';
    cfg.surffile        = 'surface_pial_left.mat';
    cfg.surfinflated    = 'surface_inflated_left.mat';
    cfg.camlight        = 'no';
    cfg.funcolormap     = flipud(mycoldmap);
    
    ft_sourceplot(cfg, int_source);
    
    view(90,0)
    h = colorbar;
    ylabel(h,'t value')
    axis off
    axis equal;
    set(gca,'fontsize',8)
    caxis([-4 0])
   
    title(sprintf('%sms',num2str(1000*mean(source_TOIs(iTOI,:)))));
    
    end
end

%%