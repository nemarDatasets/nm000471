function PNAS_F2_BC(datadir)
% PNAS_F2_BC(datadir)
% code to produce Figure 2, panels B and C
%
% datadir ........  path to 'sharing' folder (e.g., 'C:\mydata\sharing')

close all

if contains(path,'fieldtrip')
    ft_defaults;
else
    error('FieldTrip not in path');
end

addpath(fullfile(datadir,'auxiliary'));

%% settings

settings            = [];

settings.subjects   = {'s01' 's02' 's03' 's04' 's05' 's06' 's07' 's08' 's09' 's10' 's11'};
settings.conditions = {'Remember' 'Forgot'};
settings.colors     = {'g' 'r'};
settings.TRL        = [-.5 2.6];

%% preparations

load(fullfile(datadir,'iEEG',[settings.subjects{1} '_recall_TFR_HIPP']),'freq')

timeaxis = freq.time(nearest(freq.time,settings.TRL(1)):nearest(freq.time,settings.TRL(2)));
freqaxis = freq.freq;

ubermat = zeros([numel(settings.conditions) length(freqaxis) length(timeaxis) numel(settings.subjects)]);

pow_all = cell(numel(settings.subjects),numel(settings.conditions));

%% create a condition x frequency x time x subjects matrix

for isubject = 1:numel(settings.subjects)
    
    fprintf('working on subject %02d/%02d\n',isubject,numel(settings.subjects));

    load(fullfile(datadir,'iEEG',[settings.subjects{isubject} '_recall_TFR_HIPP']),'freq')

    %--- reduce to TRL time of interest

    cfg         = [];
    cfg.latency = settings.TRL;
    freq        = ft_selectdata(cfg,freq);

    %--- prepare trialinfo tags

    tmp                 = cell2mat(freq.trialinfo);
    MemoryPerformance   = [tmp.Memory];
    clear tmp
    
    %--- mark outliers within each condition
    for icondition = 1:numel(settings.conditions)
        switch icondition
            case 1 % 'Remember'
                trialselection = ismember(MemoryPerformance,{'R'});
            case 2 % 'Forgot'
                trialselection = ismember(MemoryPerformance,{'F','invalid'});
        end
        
        pow_condition                           = freq.powspctrm(trialselection==1,:,:,:);
        thresh                                  = repmat(prctile(pow_condition,90,1),[sum(trialselection) 1 1 1]);
        pow_condition(pow_condition > thresh)   = nan;
        freq.powspctrm(trialselection==1,:,:,:) = pow_condition;
        
        clear pow_condition thresh trialselection
    end
    
    %--- normalise data via z score across all time points (per channel and frequency)
    
    powz = nan(size(freq.powspctrm));
    
    for ichan = 1:size(freq.powspctrm,2)
        for ifreq = 1:size(freq.powspctrm,3)
            d = squeeze(freq.powspctrm(:,ichan,ifreq,:));
            m = nanmean(d(:));
            s = nanstd(d(:));
            clear d
            powz(:,ichan,ifreq,:) = (freq.powspctrm(:,ichan,ifreq,:) - m)./s;
            clear m s
        end
    end
    freq.powspctrm = powz;
    clear powz 
            
    for icondition = 1:numel(settings.conditions)
        
        switch icondition
            case 1 % 'Remember'
                trialselection = ismember(MemoryPerformance,{'R'});
            case 2 % 'Forgot'
                trialselection = ismember(MemoryPerformance,{'F','invalid'});
        end
        
        cfg         = [];
        cfg.trials  = find(trialselection);
        tmp         = ft_selectdata(cfg,freq);
        
        clear trialselection
        
        % store power values for later single-trial analysis
        pow_all{isubject,icondition} = tmp.powspctrm;
        
        % average across channels and trials
        ubermat(icondition,:,:,isubject) = squeeze(nanmean(nanmean(tmp.powspctrm,2),1));
        
        clear tmp
        
    end % of icondition
    clear freq MemoryPerformance
end

%% stats preparation

c_stats = cell(1,numel(settings.conditions));

for icondition=1:numel(settings.conditions)
    
    c_stats{icondition}             = [];
    c_stats{icondition}.label       = {'HIPP'};
    c_stats{icondition}.freq        = freqaxis;
    c_stats{icondition}.time        = timeaxis;
    c_stats{icondition}.dimord      = 'subj_chan_freq_time';
    c_stats{icondition}.powspctrm   = permute(ubermat(icondition,:,:,:),[4 1 2 3]);
end

%%  run stats

cfg                     = [];
cfg.latency             = settings.TRL;
cfg.frequency           = [freqaxis(1) freqaxis(end)];
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; 
cfg.correctm            = 'cluster';
cfg.alpha               = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha';
cfg.neighbours          = [];
cfg.numrandomization    = 1000;
cfg.clusterstatistic    = 'maxsum';
cfg.clusteralpha        = cfg.alpha;
cfg.clustertail         = cfg.tail;

% set up design matrix
design      = [];
design(1,:) = repmat(1:numel(settings.subjects),[1 2]);
design(2,:) = [ones(1,numel(settings.subjects)) 2*ones(1,numel(settings.subjects))];

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

Fieldtripstats = ft_freqstatistics(cfg, c_stats{:});

%% plot significant cluster (Figure 2B, left)

stats_time  = nearest(timeaxis,cfg.latency(1)):nearest(timeaxis,cfg.latency(2));
stats_freq  = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

sigmap                          = zeros(numel(freqaxis),numel(timeaxis));
sigval                          = squeeze(Fieldtripstats.mask .* Fieldtripstats.stat);
sigmap(stats_freq,stats_time)   = sigval;

figure;
uimagesc(timeaxis,freqaxis,sigmap);
axis xy
xlabel('time (s)')
xlim(settings.TRL)
ylabel('frequency (Hz)')
vline(0,'w')
title(sprintf('"%s" vs. "%s"',settings.conditions{:}));
hcb = colorbar('Location','EastOutside');
caxis([-7 7])
title(hcb,'t')

hotcold = [flipud(fliplr(colormap('hot'))); colormap('hot')];

colormap(hotcold)
set(gca,'FontName','Helvetica','FontSize',8,'TickDir','out','Box','off')

%% timecourse - average across frequencies and plot

this_FOIs = [80 90; 8 12];

for ifreq = 1:size(this_FOIs,1)
    
    foiidx    = nearest(freqaxis, this_FOIs(ifreq,1)):nearest(freqaxis, this_FOIs(ifreq,2));
    c_val_pow = nan(numel(settings.conditions),numel(timeaxis));
    
    for icondition=1:numel(settings.conditions)
        c_val_pow(icondition,:) = squeeze(nanmean(nanmean(ubermat(icondition,foiidx,:,:),2),4));
    end
    
    within_error = squeeze(nanstd(diff(nanmean(ubermat(:,foiidx,:,:),2),1,1),0,4))' ./ sqrt(size(ubermat,4));
    
    figure
    l = [];
    for icondition=1:numel(settings.conditions)
        boundedline(timeaxis,c_val_pow(icondition,:),within_error,settings.colors{icondition},'alpha')
        l(icondition) = plot(timeaxis,c_val_pow(icondition,:),settings.colors{icondition},'linewidth',2);
    end
    
    axis tight
    vline(0)
    xlabel('time (s)')
    ylabel('z power')
    title(sprintf('%d-%d Hz',this_FOIs(ifreq,:)))
    legend(l,settings.conditions);
    set(gca,'FontSize',8,'FontName','Helvetica','TickDir','out','Box','off')
    
end

%% 