function [data_wo_artifacts,IED_art_table] = rejectartifact_from_quon2022(IEDs_Quon, data,cfg_all_trials)

%% LOAD DATA
trl = cfg_all_trials.trl;
% remember the original trial definition, this may include additional columns
trlold = trl;
size_data = size(data.powspctrm); %trl, channel,freqs, time
begsample = IEDs_Quon.spikeStart_fullsampling;
endsample = IEDs_Quon.spikeEnd_fullsampling;
all_elecs = IEDs_Quon.channels;
art_seconds = []; %in what trial did it happen?
art_counter = 1;
flag_rejection = {'affected_electrodes'}; %affected_electrodes all
elec_list = [];
%% FIND ARTIFATCS AND CONVERT TO NAN
for art = 1:numel(begsample)
    %DEFINE CURRENT TRIAL
    art_trial = find(begsample(art) > trl(:,1) & begsample(art) < trl(:,2));
    if ~isempty(art_trial) & ismember(art_trial, data.trialinfo(:,8))
        %DEFINE TRIAL IV (IN FS)
        rejecttrial = trl(art_trial,1):trl(art_trial,2);
        %DEFINE TRIAL IV (IN SEC)
        trl_fromsample_tosec = linspace(min(data.time), max(data.time),length(rejecttrial));
    
        %FIND ARTIFACTUAL TIME
        time = rejecttrial>=begsample(art)  & rejecttrial<=endsample(art);
        time_sec = trl_fromsample_tosec(time); %for visual inspection
        
        %ASSIGN VALUES
        art_seconds = [art_seconds;[time_sec(1) time_sec(end) art_trial]]; 
        time_sec_iv = data.time >= art_seconds(art_counter,1) & data.time <= art_seconds(art_counter,2);
        elec_list = [elec_list; all_elecs(art)]; %#ok<*AGROW>
        
        %Remove IEDs only from electrodes showing IEDs
        current_elecs = all_elecs(art); %get names
        current_elecs = split(current_elecs, '_');
        current_elecs =  current_elecs(~cellfun('isempty',current_elecs)); %Remove empty cells
        elecs_w_IEDs = contains(data.label,current_elecs);

        %We have less trials than detected: i.e., when the spike detector
        %is used, it encompasses the whole duration of the experiment but
        %here we have less trials. Therefore, we need to index it
        %correctly.

        idx_trl_updated = find(art_trial == data.trialinfo(:,8));
        switch cell2mat(flag_rejection)
            case {'all'}
                data.powspctrm(idx_trl_updated,:,:,time_sec_iv) = nan;
            case {'affected_electrodes'}
                data.powspctrm(idx_trl_updated,elecs_w_IEDs,:,time_sec_iv) = nan;
        end
        %test functionality
        %test = squeeze(data.powspctrm(1,2,:,:));
        %test2 = squeeze(data.powspctrm(1,3,:,:));
        %test temporality
        %test3 = squeeze(data.powspctrm(78,13,:,:));
        %test4 = squeeze(data.powspctrm(78,12,:,:));
        art_counter = art_counter+1;
    end %art found
end %art
IED_art_table = array2table(art_seconds, 'VariableNames',{'Start(s)', 'End(s)', 'Trial'});
IED_art_table.elecs = elec_list;
%% PLOTTING RESULTS
figure
title('individual trials aligned to time-zero')
xlabel('time (s)');
ylabel('trial number');
hold on
% COLOR = 'grcmykgrcmykgrcmykgrcmyk';
for i=1:size(trl,1)
    %This is needed becuase some trials won't have artifacts
    x = [trl_fromsample_tosec(1) trl_fromsample_tosec(end)];
    y = [i i];
    plot(x, y, 'b')
    %Find if the current trial has spikes
    idx_current_trl = find(art_seconds(:,3) == i); %give idx of trials with spikes

    if ~isempty(idx_current_trl) %Do something to the line if some spikes present (i = art(:,3))
        for j=1:numel(idx_current_trl)
            x = trl_fromsample_tosec;
            y = ones(1,size(trl_fromsample_tosec,2))*i;
            %Find where in the trial (e.,g -2-10) the artifact is happening
            %and transform
            [~,idx_start]=min(abs(trl_fromsample_tosec-art_seconds(idx_current_trl(j),1))); val_start = trl_fromsample_tosec(idx_start);
            [~,idx_end]=min(abs(trl_fromsample_tosec-art_seconds(idx_current_trl(j),2))); val_end = trl_fromsample_tosec(idx_end);
            y(idx_start:idx_end) = i+0.4;
            %https://www.mathworks.com/matlabcentral/answers/375710-find-nearest-value-to-specific-number
           
            plot(x, y, 'b'); 
        end 
    end
end
ylim([0,size(trl,1)+1]);
sp = sprintf('Bumps are IEDs (n = %s)',num2str(height(art_seconds)));
subtitle(sp)
ax = gca; 
ax.FontSize = 16;
data_wo_artifacts = data;
end

% [~,~,idx_start] = unique(abs(trl_fromsample_tosec-art_seconds(idx_current_trl(j),1)),'stable');
% y(1,find(trl_fromsample_tosec == art_seconds(idx_current_trl(j),1)):find(trl_fromsample_tosec == art_seconds(idx_current_trl(j),2))) = nan;
% interp1(art_seconds(idx_current_trl(j),1),trl_fromsample_tosec, 'nearest')
% min(abs(trl_fromsample_tosec-art_seconds(idx_current_trl(j),1)))
