function [binned_data, bin_names] = events_to_bins(cfg,events, ivs)
%% NOTES
%events = spikes, ripples, etc
%ivs = trials, conditions, etc

%% PARAMETERS
bin_size = cfg.bin_size;
bin_time = cfg.bin_time;
bin_number = cfg.bin_number;
event_iv_fs = cfg.event_iv_fs;

%% PREPROCESSING
%Adjust ivs and events based on the desired bin output
switch bin_time
    case 's'
        events(:,1:2) = events(:,1:2)./event_iv_fs;
        ivs = ivs./event_iv_fs;
end
%Preallocate bin data
binned_data = nan(length(ivs),bin_number);

%Find duration of each event
ivs_durations = floor(ivs(:,2) - ivs(:,1));

%% MAIN FUNCTION

for iv = 1:height(ivs)
    %Subtract the start of all events from the start of all intervals
    curr_diff_ev_iv = floor(events(:,1)-ivs(iv,1));
    for bin = 1:bin_number
        %Check whether the duration of the IV is greater than the bin size
        if ivs_durations(iv) > bin_size*bin 
            events_meeting_criteria = curr_diff_ev_iv > 0 & ... %Ensure event starts after iv 
            curr_diff_ev_iv <=bin_size*bin & ...  %Ensure less than bin top
            curr_diff_ev_iv >=bin_size*(bin-1); %Ensure is more than the previous bin
            %Obtain total events
            binned_data(iv,bin) = sum(events(events_meeting_criteria,3));
            %binned_data(iv,bin) = sum(events_meeting_criteria); 
        end    
    end
end

%% FORMATTING
%Report num spikes/bin in desired format
switch bin_time
    case 's'
        binned_data = binned_data./bin_size;
end

bin_names = cell(1,length(bin_number));
for col = 1:bin_number
    bin_names(1,col) = {num2str(bin_size*col)};
end
end

%% VAULT
%To make it comparable
% bins = max(trial_data_sub.evDurAdj(2:end-1)); %Having the long baseline doesn't really make sense
% bins = floor(bins/bin_time);
% binned_rates = nan(length(event_names),bins);
% 
% for ev = 1:length(event_names)
%     events_distance_from_current_event = round((IEDs_fs_michael(:,1)-events_start_end(ev,1))/fs_events_unique);
%     for bin = 1:bins
%         %FIND AND MEET CRITERIA
%         if trial_data_sub.evDur(ev) > bin_time*bin %Ensure is in the duration of the event
%             events_meeting_criteria = events_distance_from_current_event > 0 & ... %Ensure after event
%             events_distance_from_current_event <=bin_time*bin & ...  %Ensure less than bin top
%             events_distance_from_current_event >=bin_time*(bin-1); %Ensure is more than the previous bin
%             binned_rates(ev,bin) = sum(events_meeting_criteria); 
%         end    
%     end
% end


%transforms rate to seconds
% binned_rates = binned_rates./bin_time;
% column_names = cell(1,length(bins));
% for col = 1:bins
%     column_names(1,col) = {num2str(bin_time*col)};
% end

%% PLOT RATES IN ARBITRARY TIME: 30 SECONDS
%1.) FIND IVs FOR FIRST X SECONDS + COUNT SPIKES

%Method:  %Event has to be positive (i.e., we are subtracting start of spike - event... therefore, spikes after event 
          %should yield a positive value. Additionally, we want them to be around the threshold
% arbitrary_rate_for_songs = [];
% for ev = 1:length(event_names)
%     events_distance_from_current_event = round((IEDs_fs_michael(:,1)-events_start_end(ev,1))/fs_events_unique);
%     events_meeting_criteria = events_distance_from_current_event >= 0 & events_distance_from_current_event <=rate_of_interest;
%     arbitrary_rate_for_songs = [arbitrary_rate_for_songs ; sum(events_meeting_criteria)]; %#ok<*AGROW>
% end
% %transforms rate to seconds
% arbitrary_rate_for_songs = arbitrary_rate_for_songs/rate_of_interest;