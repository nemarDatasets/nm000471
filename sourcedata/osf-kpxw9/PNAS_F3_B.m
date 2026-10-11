function PNAS_F3_B(datadir)
% PNAS_F3_B(datadir)
% code to produce Figure F3B
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
TOI_recall      = [-0.5 2.6];

%% enter subject loop

timeaxis_localiser  = TOI_localiser(1):1/200:TOI_localiser(2);
timeaxis_recall     = TOI_recall(1):1/200:TOI_recall(2);
classification      = cell(1,numel(conditionnames));

for isubject = 1:numel(subjects)
    
    fprintf('working on subject %02d/%02d\n',isubject,numel(subjects));
    
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
    cfg.latency     = TOI_recall;
    data_recall     = ft_selectdata(cfg,data_recall);
    
    % zscore across trials within run
    
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
        data_recall.trial(trialselection,:,:) = zscore(data_recall.trial(trialselection,:,:));
    end
    
    % zscore across all trials
    
    data_localiser.trial = zscore(data_localiser.trial);
    data_recall.trial    = zscore(data_recall.trial);
    
    % localiser time x recall time classification
    
    cfg                 = [];
    cfg.metric          = 'dval';
    cfg.classifier      = 'lda';
    
    for icondition = 1:numel(conditionnames)
        
        switch icondition
            case 1 % Remember trials
                classification{icondition}(:,:,:,isubject) = mv_classify_timextime(cfg, data_localiser.trial, category_localiser, data_recall.trial(trials_R==1,:,:), category_recall(trials_R==1));
            case 2 % Forgot trials
                classification{icondition}(:,:,:,isubject) = mv_classify_timextime(cfg, data_localiser.trial, category_localiser, data_recall.trial(trials_R==2,:,:), category_recall(trials_R==2));
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

for icondition = 1:numel(conditionnames)
    switch icondition
        case 1 % Remember
            dat_1.individual(:,1,:) = OSR;
            dat_2.individual(:,1,:) = SOR;
        case 2 % Forgot
            dat_1.individual(:,1,:) = OSF;
            dat_2.individual(:,1,:) = SOF;
    end
    
    design                                          = zeros(2,2*numel(subjects));
    design(1,:)                                     = [1:numel(subjects), 1:numel(subjects)];
    design(2,1:numel(subjects))                     = 1;
    design(2,numel(subjects)+1:2*numel(subjects))   = 2;
    
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.statistic           = 'depsamplesT';
    cfg.latency             = [0 2];
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
    
    figure(icondition);
    l = [];
    cnt = 0;
    hold on
    
    boundedline(timeaxis_recall,squeeze(mean(dat_1.individual)),squeeze(std(dat_1.individual)./sqrt(numel(subjects))),'b','alpha');
    cnt = cnt+1;
    l(cnt) = plot(timeaxis_recall,squeeze(mean(dat_1.individual)),'b','linewidth',2);
    
    boundedline(timeaxis_recall,squeeze(mean(dat_2.individual)),squeeze(std(dat_2.individual)./sqrt(numel(subjects))),'r','alpha');
    cnt = cnt+1;
    l(cnt) = plot(timeaxis_recall,squeeze(mean(dat_2.individual)),'r','linewidth',2);
    
    axis tight
    ylim([-.01 .01])
    vline(0)
    hline(0)
        
    % highlight significant time windows
    nanline = nan(1,numel(timeaxis_recall));
    stattime = nearest(timeaxis_recall,cfg.latency(1)):nearest(timeaxis_recall,cfg.latency(2));
    nanline(stattime(stats.mask==1)) = 0;
    cnt= cnt+1;
    hold on
    l(cnt) = plot(timeaxis_recall,nanline,'m.-','linewidth',3);
    
    legend(l,{'OS', 'SO', 'OS vs. SO'})
    xlabel('time (s)')
    ylabel('evidence')
    set(gca,'tickdir','out','box','off','FontSize',8,'FontName','Helvetica')
    title(conditionnames{icondition})
end

%%