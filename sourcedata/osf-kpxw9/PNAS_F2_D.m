function PNAS_F2_D(datadir)
% F2_D(datadir)
% code to produce Figure 2D
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
settings.TRL        = [-.5 7];
settings.FOI        = [80 90];

settings.conditions = {'Remember' 'Forgot'};

%% preparations

pow_all     = cell(numel(settings.subjects),1);
trlinfo_all = cell(numel(settings.subjects),1);

load(fullfile(datadir,'iEEG',[settings.subjects{1} '_recall_TFR_HIPP']),'freq')

timeaxis = freq.time(nearest(freq.time,settings.TRL(1)):nearest(freq.time,settings.TRL(2)));

%% create a condition x frequency x time x subjects matrix

for isubject = 1:numel(settings.subjects)
    
    fprintf('working on subject %02d/%02d\n',isubject,numel(settings.subjects));
    
    load(fullfile(datadir,'iEEG',[settings.subjects{isubject} '_recall_TFR_HIPP']),'freq')
    
    %--- reduce to TRL time and FOI frequency of interest
    cfg             = [];
    cfg.latency     = settings.TRL;
    cfg.frequency   = settings.FOI;
    freq            = ft_selectdata(cfg,freq);
    
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
    
    %--- get trial-specific power values for "Remember"
    
    trialselection = ismember(MemoryPerformance,{'R'});
    
    clear MemoryPerformance
    
    cfg         = [];
    cfg.trials  = find(trialselection);
    freq        = ft_selectdata(cfg,freq);
    
    clear trialselection
    
    %--- store power values for later single-trial analysis
    
    pow_all{isubject}       = squeeze(nanmean(nanmean(freq.powspctrm,3),2));
    trlinfo_all{isubject}   = freq.trialinfo;
    
    clear freq
end


%% look at individual trials vis a vis RTs

this_dat_all = [];
this_RT_all  = [];

for isubject = 1:numel(settings.subjects)
    
    this_trlinfo = cell2mat(trlinfo_all{isubject});
    RTs          = [this_trlinfo.RT];
    
    sel          = RTs >= settings.TRL(1) & RTs <= settings.TRL(2);
    RTs          = RTs(sel);
    this_dat     = squeeze(pow_all{isubject}(sel,:));
    
    for itrl = 1:size(this_dat,1)
        RTidx = nearest(timeaxis,RTs(itrl));
        this_dat(itrl,RTidx+1:end) = nan;
    end
        
    sel = ones(size(this_dat,1),1);
    for itrl = 1:numel(sel)
        sel(itrl) = not(any(isnan(this_dat(itrl,1:nearest(timeaxis,RTs(itrl))))));
    end
    
    this_dat    = this_dat(sel==1,:);
    RTs         = RTs(sel==1);
    
    this_dat_all = [this_dat_all;this_dat];
    this_RT_all  = [this_RT_all;RTs'];
    
end

[RTs_sorted,sort_idx] = sort(this_RT_all);

figure;
imagesc(timeaxis,[],this_dat_all(sort_idx,:));
axis xy;
colormap(jet)
hold on
plot(RTs_sorted,1:numel(RTs_sorted),'w-','linewidth',2)
vline(0)
xlabel('time (s)')
ylabel('trial #')
set(gca,'tickdir','out','box','off','FontName','Helvetica','FontSize',8)
caxis([-2 2])
colorbar

%%