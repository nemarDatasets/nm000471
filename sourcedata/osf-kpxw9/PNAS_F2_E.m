function PNAS_F2_E(datadir)
% PNAS_F2_E(datadir)
% code to produce Figure 2, panel E
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
settings.TRL        = [-1 10];
settings.FOI        = [80 90];
settings.pre        = 1.5;
settings.post       = .5;

timeaxis            = -settings.pre : .01 : settings.post;

ubermat = zeros([numel(settings.conditions) length(timeaxis) numel(settings.subjects)]);

%% create a condition x frequency x time x subjects matrix

for isubject = 1:numel(settings.subjects)
    
    fprintf('working on subject %02d/%02d\n',isubject,numel(settings.subjects));
    
    load(fullfile(datadir,'iEEG',[settings.subjects{isubject} '_recall_TFR_HIPP']),'freq')
    
    %--- reduce to TRL time and FOI frequency of interest    
    cfg           = [];
    cfg.latency   = settings.TRL;
    cfg.frequency = settings.FOI;
    freq          = ft_selectdata(cfg,freq);
    
    %--- prepare trialinfo tags
    tmp                 = cell2mat(freq.trialinfo);
    MemoryPerformance   = [tmp.Memory];
    RTs                 = [tmp.RT];
    
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
        
        % collapse across channels
        dat         = squeeze(nanmean(tmp.powspctrm,2));
        cond_RTs    = RTs(trialselection==1);
        
        clear trialselection
        
        %-- cut out time
        
        dat_tw = nan([size(dat,1),size(dat,2),length(timeaxis)]);
        
        for itrl = 1:size(dat,1)
            
            rt = cond_RTs(itrl);
            if not(isfinite(rt))
                rt=-1;
            end
            rt = round(rt*100)/100;
            if rt < 0 % invalid trial, skip
                disp(['trl:' num2str(itrl) ' rt:' num2str(rt)])
                
            elseif (rt - settings.pre) <= 0 && (rt + settings.post) < freq.time(end) % RT-settings.pre out of bounds, RT+settings.post within bounds
                disp(['trl:' num2str(itrl) ' rt:' num2str(rt)])
                
                idx = nearest(freq.time,0):nearest(freq.time,rt+settings.post);
                
                dat_tw(itrl,:,(length(timeaxis)-length(idx)+1):length(timeaxis)) = squeeze(dat(itrl,:,idx));
                
            elseif (rt - settings.pre) > 0 && (rt + settings.post) > freq.time(end) % RT-settings.pre within bounds, RT+settings.post out of bounds
                disp(['trl:' num2str(itrl) ' rt:' num2str(rt)])
                
                idx = nearest(freq.time,rt-settings.pre):nearest(freq.time,rt+settings.post);
                
                dat_tw(itrl,:,1:length(idx)) = squeeze(dat(itrl,:,idx));
                
                
            else % RT-settings.pre : RT+settings.post are within available time range
                
                idx = nearest(freq.time,rt-settings.pre):nearest(freq.time,rt+settings.post);
                
                dat_tw(itrl,:,:) = squeeze(dat(itrl,:,idx));
            end
        end % of itrl
        
        % average across channels and trials
        ubermat(icondition,:,isubject) = squeeze(nanmean(nanmean(dat_tw,2),1));
        
        clear tmp
        
    end % of icondition
    clear freq MemoryPerformance
end

%% timecourse

c_val_pow = nan(numel(settings.conditions),numel(timeaxis));

for icondition=1:numel(settings.conditions)
    c_val_pow(icondition,:) = squeeze(nanmean(ubermat(icondition,:,:),3));
end

within_error = squeeze(nanstd(diff(ubermat,1,1),0,3))' ./ sqrt(size(ubermat,3));

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
title(sprintf('%d-%d Hz',settings.FOI))
legend(l,settings.conditions);
set(gca,'FontSize',8,'FontName','Helvetica','TickDir','out','Box','off')

%%