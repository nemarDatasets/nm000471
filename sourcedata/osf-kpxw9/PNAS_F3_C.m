function PNAS_F3_C(datadir)
% PNAS_F3_C(datadir)
% code to produce Figure F3C
%
% datadir ... path to 'sharing' folder (e.g., 'C:\mydata\sharing')

close all

if contains(path,'fieldtrip')
    ft_defaults;
else
    error('FieldTrip not in path');
end

addpath(fullfile(datadir,'auxiliary'));

%% settings

subjects        = {'s01' 's02' 's03' 's04' 's05' 's06' 's07' 's08' 's09' 's10' 's11'};
conditionnames  = {'Remember' 'Forgot'};
TOI_localiser   = [0 1];
TOI_peaksearch  = [0 2.6];
TOI_around_peak = [-.5 .5];
TOI_TFR_recall  = [-.5 2.6];
TOI_EEG_recall  = TOI_peaksearch + TOI_around_peak; % for classification, needs to be peak search range plus/minus TOI_around_peak

%% enter subject loop

timeaxis_localiser  = TOI_localiser(1):1/200:TOI_localiser(2);
timeaxis_recall     = -.5:1/200:.5;
classification      = cell(1,numel(conditionnames));

for isubject = 1:numel(subjects)
    
    fprintf('working on subject %02d/%02d\n',isubject,numel(subjects));
    
    %% gamma peaks
    
    load(fullfile(datadir,'iEEG',[subjects{isubject} '_recall_TFR_HIPP']),'freq')
    
    %--- reduce to TRL time of interest
    cfg             = [];
    cfg.latency     = TOI_TFR_recall;
    cfg.frequency   = [80 90];
    freq            = ft_selectdata(cfg,freq);
    
    %--- prepare trialinfo tags
    tmp                 = cell2mat(freq.trialinfo);
    MemoryPerformance   = [tmp.Memory];
    clear tmp
    
    %--- mark outliers within each condition
    for icondition = 1:numel(conditionnames)
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
    
    pow = squeeze(nanmean(nanmean(freq.powspctrm,2),3));
    
    %% search single-trial power peaks
    
    peak_window = nearest(freq.time,TOI_peaksearch(1)):nearest(freq.time,TOI_peaksearch(2));
    
    peaks = nan(size(pow,1),1);
    
    for itrial=1:size(pow,1)
        d = pow(itrial,peak_window);
        [pks,locs] = findpeaks(pow(itrial,peak_window));
        [~,maxpk] = max(pks);
        ix = locs(maxpk);
        if not(isempty(ix))
            peaks(itrial) = freq.time(peak_window(ix));
        end
    end
    
    %% classification
    
    data_localiser  = load(fullfile(datadir,'iEEG',[subjects{isubject} '_localiser_EEG_extraHIPP.mat']));
    data_localiser  = data_localiser.eeg;
    
    data_recall     = load(fullfile(datadir,'iEEG',[subjects{isubject} '_recall_EEG_extraHIPP.mat']));
    data_recall     = data_recall.eeg;
    
    % bring data into matrix format
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.removemean  = 'no';
    data_localiser  = ft_timelockanalysis(cfg, data_localiser);
    data_recall     = ft_timelockanalysis(cfg, data_recall);
    
    % get trial info for classification
    % localiser: object cue = 1,  scene cue = 2
    
    trlinfo_localiser  = cell2mat(data_localiser.trialinfo);
    category_localiser = [trlinfo_localiser.StimType];
    
    % recall
    trlinfo_recall     = cell2mat(data_recall.trialinfo);
    
    % recall: 'Remember'  = 1, 'Forgot' = 2
    trials_R                                = double(ismember([trlinfo_recall.Memory],'R'));
    trials_R(not(trials_R))                 = 2;
    
    % recall: object cue = 1,  scene cue = 2
    category_recall                         = double([trlinfo_recall.ObjectProbe]);
    category_recall(not(category_recall))   = 2;
    
    % apply 100 ms smoothing window to localiser and recall data
    data_localiser.trial = smoothdata(data_localiser.trial,3,'movmean',.100/(1/200));
    data_recall.trial    = smoothdata(data_recall.trial,3,'movmean',.100/(1/200));
    
    % reduce to TOIs
    cfg             = [];
    cfg.latency     = TOI_localiser;
    data_localiser  = ft_selectdata(cfg,data_localiser);
    
    cfg             = [];
    cfg.latency     = TOI_EEG_recall;
    data_recall     = ft_selectdata(cfg,data_recall);
    
    %% realign retrieval data based on gamma peaks
    
    fs = mode(diff(data_recall.time));
    
    ret_realigned = nan(size(data_recall.trial,1),numel(data_recall.label),round(diff(TOI_around_peak)/fs)+1);
    
    valid_trials = nan(size(data_recall.trial,1),1);
    for itrial=1:size(data_recall.trial,1)
        valid_trials(itrial) = isfinite(peaks(itrial));
        if valid_trials(itrial)
            idx = nearest(data_recall.time,peaks(itrial))+round(TOI_around_peak(1)/fs):nearest(data_recall.time,peaks(itrial))+round(TOI_around_peak(2)/fs);
            ret_realigned(itrial,:,:) = data_recall.trial(itrial,:,idx);
        end
    end
    
    ret_realigned           = ret_realigned(valid_trials==1,:,:);
    trials_R                = trials_R(valid_trials==1);
    category_recall         = category_recall(valid_trials==1);
    trlinfo_recall          = trlinfo_recall(valid_trials==1);
    
    %% zscore across trials within run
    
    localiser_run_info   = [trlinfo_localiser.RunNumber];
    localiser_runs_total = unique(localiser_run_info);
    
    for irun=1:numel(localiser_runs_total)
        trialselection = localiser_run_info == localiser_runs_total(irun);
        data_localiser.trial(trialselection,:,:) = zscore(data_localiser.trial(trialselection,:,:));
    end
    
    recall_run_info   = [trlinfo_recall.RunNumber];
    recall_runs_total = unique(recall_run_info);
    
    for irun=1:numel(recall_runs_total)
        trialselection = recall_run_info == recall_runs_total(irun);
        ret_realigned(trialselection,:,:) = zscore(ret_realigned(trialselection,:,:));
    end
    
    % zscore across all trials
    data_localiser.trial = zscore(data_localiser.trial);
    ret_realigned    = zscore(ret_realigned);
    
    % localiser time x recall time classification
    cfg                 = [];
    cfg.metric          = 'dval';
    cfg.classifier      = 'lda';
    
    for icondition = 1:numel(conditionnames)
        
        switch icondition
            case 1 % Remember trials
                classification{icondition}(:,:,:,isubject) = mv_classify_timextime(cfg, data_localiser.trial, category_localiser, ret_realigned(trials_R==1,:,:), category_recall(trials_R==1));
            case 2 % Forgot trials
                classification{icondition}(:,:,:,isubject) = mv_classify_timextime(cfg, data_localiser.trial, category_localiser, ret_realigned(trials_R==2,:,:), category_recall(trials_R==2));
        end
    end
end

%% timecourses

OSR = [];
OSF = [];
SOR = [];
SOF = [];

for ii=1:numel(subjects)
    
    OSR(ii,1,:,:) = squeeze(classification{1}(nearest(timeaxis_localiser,.35),1,:,ii));
    OSF(ii,1,:,:) = squeeze(classification{2}(nearest(timeaxis_localiser,.35),1,:,ii));
    SOR(ii,1,:,:) = squeeze(classification{1}(nearest(timeaxis_localiser,.35),2,:,ii));
    SOF(ii,1,:,:) = squeeze(classification{2}(nearest(timeaxis_localiser,.35),2,:,ii));
    
end

tmp        = struct();
tmp.time   = timeaxis_recall;
tmp.label  = {'LDA'};
tmp.dimord = 'subj_chan_time';

dat_1 = tmp;
dat_2 = tmp;

dat_1.individual(:,1,:) = OSR;
dat_2.individual(:,1,:) = SOR;

design                                          = zeros(2,2*numel(subjects));
design(1,:)                                     = [1:numel(subjects), 1:numel(subjects)];
design(2,1:numel(subjects))                     = 1;
design(2,numel(subjects)+1:2*numel(subjects))   = 2;

cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.statistic           = 'depsamplesT';
cfg.latency             = TOI_around_peak;
cfg.clusterstatistic    = 'maxsum';
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.correcttail         = 'alpha';
cfg.tail                = 0;
cfg.clustertail         = cfg.tail;
cfg.numrandomization    = 1000;
cfg.neighbours          = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;

stats                   = ft_timelockstatistics(cfg, dat_1,dat_2);

figure;
l = [];
cnt=0;
hold on
boundedline(timeaxis_recall,squeeze(mean(dat_1.individual)),squeeze(std(dat_1.individual)./sqrt(numel(subjects))),'b','alpha');
cnt= cnt+1;
l(cnt) = plot(timeaxis_recall,squeeze(mean(dat_1.individual)),'b','linewidth',2);

boundedline(timeaxis_recall,squeeze(mean(dat_2.individual)),squeeze(std(dat_2.individual)./sqrt(numel(subjects))),'r','alpha');
cnt= cnt+1;
l(cnt) = plot(timeaxis_recall,squeeze(mean(dat_2.individual)),'r','linewidth',2);

axis tight
ylim([-.005 .005])
vline(0)
hline(0)

% highlight significant time windows
nanline  = nan(1,numel(timeaxis_recall));
stattime = nearest(timeaxis_recall,cfg.latency(1)):nearest(timeaxis_recall,cfg.latency(2));
nanline(stattime(stats.mask==1)) = 0;
cnt= cnt+1;
hold on
l(cnt) = plot(timeaxis_recall,nanline,'m.-','linewidth',3);

legend(l,{'OS', 'SO', 'OS vs. SO'})
xlabel('time (s)')
ylabel('evidence')
set(gca,'tickdir','out','box','off','FontSize',8,'FontName','Helvetica')
title('Remember')

%%