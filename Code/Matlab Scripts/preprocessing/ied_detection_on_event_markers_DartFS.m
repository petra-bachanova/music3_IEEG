function cfg_out_IEDs_with_events = ied_detection_on_event_markers_DartFS(cfg_IEDs_events, subjID)

cwd = pwd;

%% NOTES:

%CONSTANT VALUES IN THE EXPERIMENT
% self.BASELINE_DUR = 60        # 60
% self.COUNTDOWN_TO_SART = 5    # 5
% self.TEXT_HEIGHT = 100        # 100
% 
% self.SART_DIGITS = 4        # 4
% self.SART_DIGIT_DUR = 0.5   # 1.0
% self.SART_MASK_DUR = 0.2    # 0.2
% self.SART_DUR = 60          # 60
% self.SART_NUMBER_OF_CIRCLES = int(self.SART_DUR / (self.SART_DIGIT_DUR + self.SART_MASK_DUR))

%COMPUTING EVENT DURATIONS
%Trial Structure: 1-minute baseline || Music || Pause || SART || 1-minute baseline 

%Casey based:
%1.) Estimate 1-minute baseline for every song
%2.) Estimate SART to be 1 second
%3.) Estimate Pause to be anything remaining

%Log-based:
%1.) Compute duration of pause and sart based on log file
%fil accordingly


%% CREATE FOLDER
folder = strcat('/Volumes/Ecog/music3_IEEG/Results/ied_rates/', subjID);
if ~exist(folder, 'dir')
       mkdir(folder)
end
cd(folder)


%% PARSE ARGUMENTS
IEDs_Quon = cfg_IEDs_events.IEDs_Quon.IEDs; 
trial_data = cfg_IEDs_events.trial_data;
SART_data = cfg_IEDs_events.SART_data;
subjID = cfg_IEDs_events.subjID;
MusicSession = cfg_IEDs_events.MusicSession;
hdr = cfg_IEDs_events.hdr;
rate_of_interest = cfg_IEDs_events.rate_of_interest;
bin_time = cfg_IEDs_events.bin_time;
ev_times_data = readtable(['/Volumes/Ecog/music3_IEEG/Data/Log_files/' subjID '_event_times.xlsx']);

%PREPROCESS ARGS
last_time_eegdata = hdr.nSamples/hdr.Fs; % EXPERIMENT DURATION
trial_data_sub = trial_data(contains(trial_data.subject,subjID, "IgnoreCase",true) & trial_data.MusicSession == MusicSession,:);
behavioral_data_sub = SART_data(contains(SART_data.subject,subjID) & SART_data.MusicSession == MusicSession,:);
behavioral_data_evs = unique(behavioral_data_sub.EV,'stable');
fs_events = trial_data_sub.sfreq; 
fs_events_unique = unique(fs_events);
last_time_eegdata_michael = round(last_time_eegdata*fs_events_unique);


%% LOAD DATA
events_start_end = [trial_data_sub.startSample trial_data_sub.endSample];
song_durations = diff(events_start_end,1,2);
event_names = trial_data_sub.EV;
song_names = event_names; %In michael's format
event_names = replace(event_names, '_', '-'); %for vis

if numel(unique(fs_events)) > 1 %This can be easily fixed my multiplying each event by its fs
    warning('Fs is not consistent... code does not support')
end


%% ADD EVENTS

%TRIAL DATA
%Expand table
newlength = height(trial_data_sub) + (height(events_start_end))*2 + 2; %SART + Baseline per each event + long baseline + pad
trial_data_sub(newlength,:) = trial_data_sub(1,:); %append event in the last row
trial_data_sub(end,:)=[]; %delete padd (allowed to match data)

%Fill in repeated data
trial_data_sub.subject(height(events_start_end)+1:end) = trial_data_sub.subject(1);
trial_data_sub.MusicSession(height(events_start_end)+1:end) = trial_data_sub.MusicSession(1);

num_events = numel(event_names);
adding_counter = 1+ num_events;

if num_events > numel(behavioral_data_evs)
    warning('The log file doesnt have enough data to determine the SART')
end

% 1- MINUTE BASELINES & SART
if ~any(contains(event_names, 'baseline'))
    for ev = 1:num_events
        %BASELINE (60 seconds)
        event_names(end+1) = strcat(event_names(ev),'-baseline'); %Add name
        %Append event time at the end of matrix
        events_start_end = [events_start_end; ([events_start_end(ev,1) - fs_events(ev)*60 events_start_end(ev,1)-1]) ]; 
        %Append to trial data
        trial_data_sub.EV(adding_counter) = event_names(adding_counter);
        trial_data_sub.startSample(adding_counter) = events_start_end(adding_counter,1);
        trial_data_sub.endSample(adding_counter) = events_start_end(adding_counter,2);
        trial_data_sub.sfreq(adding_counter) = fs_events(ev);
        trial_data_sub.logFEv(adding_counter) = adding_counter; %Could be better
        trial_data_sub.evDur(adding_counter) = (trial_data_sub.endSample(adding_counter) - trial_data_sub.startSample(adding_counter))/trial_data_sub.sfreq(adding_counter);
        adding_counter = adding_counter+1;
   


        %DEFINE DURATION OF SART AND BASELINE
        event_names(end+1) = strcat(event_names(ev),'-SART');
        behavioral_events_computation = 'michael';

        switch behavioral_events_computation
            case 'log_file'
                %SONG DURATION
                sound_end = find(contains(ev_times_data.EV, [cell2mat(song_names(ev)), '.wav -- Stoptime']), 1,'first');
                SART = find(contains(ev_times_data.EV, 'SART'));
                idx_start_current_ev = find(contains(ev_times_data.EV, song_names(ev)), 1,'first');
                idx_end_current_ev =   find(contains(ev_times_data.EV, song_names(ev)), 1,'last');
                if idx_start_current_ev ~= idx_end_current_ev
                    %Event duration computed from 
                    duration_current_ev = str2double(ev_times_data.start_time(idx_end_current_ev)) - str2double(ev_times_data.start_time(idx_start_current_ev));
                elseif idx_start_current_ev == idx_end_current_ev
                    warning('The duration of the song cant be determined. We will use Michaels data')
                end
                %SART
                if isempty(sound_end) || ev > numel(SART)
                    warning('We dont have data to estimate the gap between the pause and start of new trial.')
                    %Question: NaN or estimation
                    events_start_end = [events_start_end; ([events_start_end(ev,2)+floor(fs_events(ev)*10) ... %take end of event
                        events_start_end(ev,2)+floor(fs_events(ev)*10)+floor(fs_events(ev)*60)])];
                    %events_start_end = [events_start_end; [NaN NaN]];

                else
                    SART = SART(ev);
                    gap_sound_pause = str2double(ev_times_data.start_time(SART)) - str2double(ev_times_data.start_time(sound_end));

                    %Find duration of SART
                    idx_start_current_SART = find(contains(behavioral_data_sub.EV,behavioral_data_evs(ev)), 1,'first');
                    idx_end_current_SART = find(contains(behavioral_data_sub.EV,behavioral_data_evs(ev)), 1,'last');
                    duration_current_SART = str2double(behavioral_data_sub.AppearanceTime(idx_end_current_SART)) - str2double(behavioral_data_sub.AppearanceTime(idx_start_current_SART));
                    events_start_end = [events_start_end; ([events_start_end(ev,2)+floor(fs_events(ev)*gap_sound_pause) ... %take end of even
                        events_start_end(ev,2)+floor(fs_events(ev)*gap_sound_pause)+floor(fs_events(ev)*duration_current_SART)])];
                end
                % sound_end = find(contains(SART_data.EV), '-- stoptime'), ev);
                % pause = find(contains(SART_data.EV), '-- PAUSE'), ev);

            case 'michael'
                %SONG DURATION
                duration_current_ev = song_durations(ev)/fs_events(ev);
                %SART: (1-minute baseline || Music || Pause (?) || SART (60) || 1-minute baseline)
                
                %We can't reliably estimate SART duration for last event
                %using Micheal's data only
                if ev<num_events
                    pause_and_SART = events_start_end(ev+1,1) - events_start_end(ev,2) ... %Subtract the start of song 2 - song 1 = Pause + SART + baseline
                        -  fs_events_unique*60; %subtract baseline
                elseif ev==num_events
                    pause_and_SART = fs_events_unique*70;
                end

                if pause_and_SART/fs_events_unique < 60
                    warning('The pause and baseline shouldnt be less than 60 seconds')
                    events_start_end = [events_start_end; ([events_start_end(ev,2)+1 ... 
                                  events_start_end(ev,2) + duration_current_ev])];

                elseif pause_and_SART/fs_events_unique > 60
                    pause_duration = pause_and_SART - 60*fs_events_unique;
                    events_start_end = [events_start_end; ([events_start_end(ev,2)+1+pause_duration ... 
                                  events_start_end(ev,2)+1+pause_and_SART])];
                end
               
                
                %2.) Estimate SART to be 60 second
                %3.) Estimate Pause to be anything remaining


        end
        
        %Append to trial data
        trial_data_sub.EV(adding_counter) = event_names(adding_counter);
        trial_data_sub.startSample(adding_counter) = events_start_end(adding_counter,1);
        trial_data_sub.endSample(adding_counter) = events_start_end(adding_counter,2);
        trial_data_sub.sfreq(adding_counter) = fs_events(ev);
        trial_data_sub.logFEv(adding_counter) = adding_counter; %Could be better
        trial_data_sub.evDur(adding_counter) = (trial_data_sub.endSample(adding_counter) - trial_data_sub.startSample(adding_counter))/trial_data_sub.sfreq(adding_counter);
        if isnan(trial_data_sub.evDur(adding_counter))
            trial_data_sub.evDur(adding_counter) = 60;
            warning('Event has been assigned to 60 seconds)')
        end
        adding_counter = adding_counter+1;
    end
end

% LONG BASELINE (COULD BE BETTER)
if ~any(contains(event_names, 'long-baseline'))
    event_names(end+1) = {strcat('long-baseline')};
    trial_data_sub.EV(end) = event_names(end);
    trial_data_sub.sfreq(end) = fs_events(1);
    trial_data_sub.logFEv(end) = adding_counter; 
    trial_data_sub = sortrows(trial_data_sub,{'startSample'},{'ascend'});
    events_start_end = [events_start_end; [0 trial_data_sub.startSample(2)-1]]; 
    trial_data_sub.startSample(1) = events_start_end(end,1);
    trial_data_sub.endSample(1) = events_start_end(end,2);
    trial_data_sub.evDur(1) = (trial_data_sub.endSample(1) - trial_data_sub.startSample(1))/trial_data_sub.sfreq(adding_counter);
end
%current_sub(:,setdiff(current_sub.Properties.VariableNames, {("startSample"), ("endSample")}))

%% SANITY CHECK: ARE THERE DIFFERENCES IN TIMING
events_start_end = [trial_data_sub.startSample trial_data_sub.endSample]; %Grab sorted time
behavioral_offset = zeros(height(events_start_end),1);
for a = 1:height(events_start_end)-1
    current_diff = events_start_end(a+1,1) - events_start_end(a,2);
    if current_diff<0
        warning('Discrepancies between Michael and log file data')
        behavioral_offset(a) = events_start_end(a,2)-events_start_end(a+1);
        events_start_end(a,2) = events_start_end(a+1);
    end
end

trial_data_sub.startSample = events_start_end(:,1);
trial_data_sub.endSample = events_start_end(:,2);
trial_data_sub.behavioral_offset = behavioral_offset;
trial_data_sub.evDurAdj = (trial_data_sub.endSample - trial_data_sub.startSample)/trial_data_sub.sfreq(1);

% MAKE START AND END OF SPIKES TO COINCIDE WITH MICHAEL'S
IEDs_2048 = [IEDs_Quon.spike_start_2048_hz IEDs_Quon.spike_end_2048_hz]; %DATA ON 2048fs
IEDs_fs_michael = convert_ivs_to_other_fs(IEDs_2048,2048, fs_events_unique);

IEDs_Quon.spike_start_michael_hz = IEDs_fs_michael(:,1);
IEDs_Quon.spike_end_michael_hz = IEDs_fs_michael(:,2);

%% ANALYSIS: SPIKES PER DETECTION
events = [IEDs_fs_michael IEDs_Quon.num_spikes];
labeled_events = label_events_in_ivs(events, events_start_end);
IEDs_Quon.part = labeled_events(:,4);
num_events = numel(event_names);

summary_table = []; %ID, total spikes
for stimuli = 1:num_events 
    summary_table(stimuli,:) = sum(IEDs_Quon.num_spikes(IEDs_Quon.part == stimuli)); 
end
%Add to table
trial_data_sub.Count = summary_table;

event_names = replace(trial_data_sub.EV, '_', '-');
summary_table = array2table(summary_table,"RowNames",event_names,'VariableNames', {'Count'});
% PREPROCESS RATES: FIND EVENT DURATIONS
event_duration_seconds = round(trial_data_sub.evDurAdj);

%%% ACTIVATE FOR GOOD PLOTS

% PLOT TOTAL QUANTITY
plot_events_totals(event_names, summary_table.Count)
title(strcat('IEDs count (Quon 2022) for subject:', subjID, ', music session:', num2str(trial_data_sub.MusicSession(1))))
saveas(gcf, ['/Volumes/Ecog/music3_IEEG/Results/ied_rates/', subjID,'_s',num2str(MusicSession),'_count.jpg'])
pause(0.3)
close

% PLOT RATES
IED_rate = plot_events_rate_seconds(summary_table.Count, event_duration_seconds, event_names);
title(strcat('IEDs count (Quon 2022) for subject:', subjID, ', music session:', num2str(trial_data_sub.MusicSession(1))))
trial_data_sub.Rate = IED_rate;
trial_data_sub.spike_baseline = repmat(trial_data_sub.Rate(1), height(trial_data_sub), 1); %Lonf baseline
saveas(gcf, ['/Volumes/Ecog/music3_IEEG/Results/ied_rates/', subjID,'_s',num2str(MusicSession),'_rate(hz).jpg'])
pause(0.3)
close


% PLOT IEDs OVER TIME: COOL GRAPH!
%plot_ieds_overtime(eegdata, trial_data_sub,IEDs_Quon)
%% SPIKE RATE IN 15S NON-OVERLAPPING BINS
cfg = [];
cfg.bin_size = 15;
cfg.bin_time = 's';
cfg.bin_number = 10;
cfg.event_iv_fs = fs_events_unique;
[binned_data, bin_names] = events_to_bins(cfg,labeled_events, events_start_end);

initial_rate_for_songs_table = array2table(binned_data,"RowNames",event_names,'VariableNames',bin_names);
max_rate = max(binned_data,[],'all');
trial_data_sub = [trial_data_sub initial_rate_for_songs_table];

% PLOT OVER DURATION OF MINIMUM EVENT
%plot_rate_in_bins(event_duration_seconds, bin_time, binned_data)

%% SPIKE RATES: FIRST 30 SECONDS
cfg = [];
cfg.bin_size = 30;
cfg.bin_time = 's';
cfg.bin_number = 1;
cfg.event_iv_fs = fs_events_unique;
[binned_data_30s, bin_names_30s] = events_to_bins(cfg,labeled_events, events_start_end);
rate_30s = array2table(binned_data_30s,"RowNames",event_names,'VariableNames', bin_names_30s);
plot_events_totals(event_names, binned_data_30s);
title(strcat('IEDs count (Quon 2022) for subject:', subjID, ', music session:', num2str(trial_data_sub.MusicSession(1)), ', First(s):', num2str(rate_of_interest)))
saveas(gca,['/Volumes/Ecog/music3_IEEG/Results/ied_rates/' subjID,'_s',num2str(MusicSession),'_30s_rate(hz).jpg'])
pause(0.3)
close
%% ADD TO TABLE
writetable(trial_data_sub, ['/Volumes/Ecog/music3_IEEG/Results/ied_rates/' subjID '_ied_rates.xlsx'])
writetable(IEDs_Quon,      ['/Volumes/Ecog/music3_IEEG/Results/ied_rates/' subjID '_ied_detection_summary.xlsx'])

%trial_data = [trial_data_sub; trial_data(~contains(trial_data.subject,subjID) & trial_data.MusicSession == MusicSession,:)];

%% SAVE
cfg_out_IEDs_with_events = struct;
cfg_out_IEDs_with_events.subjID = subjID;
cfg_out_IEDs_with_events.IED_table = IEDs_Quon;
cfg_out_IEDs_with_events.count_table = summary_table;
cfg_out_IEDs_with_events.rate_table = rate_table;
cfg_out_IEDs_with_events.arbrate_table = arbrate_table;
cfg_out_IEDs_with_events.arbrate = rate_of_interest;
cfg_out_IEDs_with_events.initial_rate_for_songs_table= initial_rate_for_songs_table;
cfg_out_IEDs_with_events.current_subject = trial_data_sub;
cfg_out_IEDs_with_events.trial_data = trial_data;

cfg_out_IEDs_with_events.part_legend = {['For the part denomination, 1 = baseline, that is, anything before the first piece of music' ...
    'the last number corresponds to ISI. That is, anything between songs']};
end

%% PLOT RATE ACROSS MIN NUMBER OF BINS
function plot_rate_in_bins(event_duration, bin_time, binned_data)
%minimum_event_duration = min(event_duration);
%minimum_complete_events = floor(minimum_event_duration/bin_time);
data_plot = binned_data(:,1:4);
figure(Position=[0 0 1500,1500])
for plt = 1:height(initial_rate_for_songs_table)
    subplot(height(initial_rate_for_songs_table),1,plt)
    plot(str2double(column_names(1:minimum_complete_events)),data_plot(plt,:),'-o','Color',COLOR.blue, 'LineWidth',2)
    ylim([0 max_rate])
    text(str2double(column_names(1:minimum_complete_events)),data_plot(plt,:),num2cell(round(data_plot(plt,:),2)),'VerticalAlignment','top','HorizontalAlignment','left')
    ax = gca;
    ax.FontSize = 16;
    set(gca,'fontname','SansSerif')
    set(gcf,'color','w')
    xticks([])
    if plt == ceil(height(initial_rate_for_songs_table)/2)
        ylabel('IED rate (Hz)')
    end
    title(event_names(plt))
end
sgtitle(strcat('IEDs count (Quon 2022) for subject:', subjID, ', music session:', num2str(trial_data_sub.MusicSession(1)), ', Binned'))
%xticklabels(str2double(column_names))
xticks(str2double(column_names))
xlabel('Time bins')
ax = gca;
ax.FontSize = 16;
set(gca,'fontname','SansSerif')
%saveas(gca,[subjID,'_s',num2str(MusicSession),'_binned_rate(hz).jpg'])
close

end
%%
function plot_ieds_overtime(eegdata, trial_data_sub,IEDs_Quon)
average_EEG = mean(eegdata.trial{1,1},1,"omitnan"); %#ok<*UNRCH>
figure(Position=[0 0 1500,1500])
%Average iEEG
subplot(2,1,1)
reduce_plot(eegdata.time{1,1},average_EEG,"Color",COLOR.black); hold on; %Plot is in seconds
for events = 1:height(trial_data_sub)
    %xline(current_sub.startSample(events)/fs_michael_unique, '--', Color=COLOR.red)
    xregion(trial_data_sub.startSample(events)/fs_michael_unique,trial_data_sub.endSample(events)/fs_michael_unique,FaceColor = COLOR.frontal3);
end
axis off; axis tight;

%SPIKE RASTER
spike_raster = zeros(numel(eegdata.label),round(eegdata.time{1,1}(end)));

for spike = 1:height(IEDs_Quon)
    current_chns = split(IEDs_Quon.channels(spike),"_");
    current_chns = current_chns(~cellfun(@isempty,current_chns));
    for idx = 1:numel(current_chns)
        spike_raster(strcmp(eegdata.label,current_chns(idx)),IEDs_Quon.spike_start_seconds(spike)) = 1;
    end
    %spike_raster(1,spikes_histo_data(spike,4)) = spikes_histo_data(spike,1);
end
subplot(2,1,2)
imagesc([0:round(eegdata.time{1,1}(end))], [1:numel(eegdata. label)], spike_raster, [0 1]);
colormap(flipud(gray(256)));
%colormap((gray(256)));
xlabel('Time (s)'); ylabel('Electrode ID')
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLength = [0 0];
end

%% PLOT EVENTS TOTALS
function plot_events_totals(event_names, count)

%Regular bar plot
X = categorical(event_names);
Y = count; 
% Sort by decreasing y value.
[sortedY_count, sortOrder_count] = sort(Y, 'descend');
sortedX_count = X(sortOrder_count);
figure(Position=[0 0 1500,1500])
barh(sortedY_count)
yticks([1:1:numel(sortedX_count)])
yticklabels(sortedX_count)
if any(mod(count(:),1) ~= 0)
    xlabel('Rate (Hz)')
else
    xlabel('Count'); 
end
ax = gca;
ax.FontSize = 16;
set(gca,'fontname','SansSerif') 
set(gcf,'color','w') %Remove ugly gray background
%saveas(gca,[subjID,'_s',num2str(MusicSession),'_total_count.jpg'])
end

%% PLOT RATES IN SECONDS
function IED_rate = plot_events_rate_seconds(ev_count, ev_duration, ev_names)

IED_rate = ev_count./ev_duration;
rate_table = array2table(IED_rate,"RowNames",ev_names,'VariableNames', {'Rate (Hz)'});

figure(Position=[0 0 1000,500])
[sortedY_rate, sortOrder_rate] = sort(IED_rate, 'descend');
sortedX_rate= ev_names(sortOrder_rate); 
barh(sortedY_rate)
yticks([1:1:numel(sortedX_rate)])
yticklabels(sortedX_rate)
xlabel('Rate (Hz)'); 
ax = gca;
ax.FontSize = 16;
set(gca,'fontname','SansSerif') 
set(gcf,'color','w') %Remove ugly white background
%saveas(gca,[subjID,'_s',num2str(MusicSession),'_rate(hz).jpg'])
end

%WEIRD: IN SOME CASES, FILE DURATION IS SHORTER THAN LAST SPIKE... CHECK LATER

% if last_time_eegdata_michael > IEDs_fs_michael(end,2)
%     last_sample = last_time_eegdata_michael/fs_events_unique;
% elseif last_time_eegdata_michael < IEDs_fs_michael(end,2)
%     last_sample = IEDs_fs_michael(end,2)/fs_events_unique;
% end
% ISI_duration = round(last_sample - sum(event_duration_seconds));
% %numel(eegdata.time{1, 1})> round((IEDs_2048(end,2)/2048)*512)
% event_duration_seconds = [event_duration_seconds ; ISI_duration];
