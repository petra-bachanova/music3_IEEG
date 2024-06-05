function labeled_events = label_events_in_ivs(events, ivs)

labeled_events = zeros(height(events),1);
for ev = 1:height(events)
    evnts_in_iv = events(ev,1) > ivs(:,1) & events(ev,2) < ivs(:,2);
    evnts_in_iv_idx = find(evnts_in_iv); %note music idx is 1 more than logFEv due to Python idx 0 and matlab 1
    if ~isempty(evnts_in_iv_idx) %) = location_idx;
        labeled_events(ev,1) = evnts_in_iv_idx; 
        %labeled_events(ev,2) = events(ev,3);
    elseif isempty(evnts_in_iv_idx) %present in the intervals
        labeled_events(ev,1) = nan;
    end
end
%Append to previous structure
labeled_events = [events labeled_events];

end
%% VAULT 
% num_events = numel(event_names); 
% spikes_per_detection = IEDs_Quon.num_spikes;
% IED_per_stimuli = zeros(height(IEDs_Quon),1); %Corresponds to number of spikes x idx of spike
% 
% for spike = 1:height(IEDs_fs_michael) 
%     %Try events that are NOT intevals
%     location_logical = IEDs_fs_michael(spike,1) > events_start_end(:,1) & IEDs_fs_michael(spike,2) < events_start_end(:,2);
%     location_idx = find(location_logical); %note music idx is 1 more than logFEv due to Python idx 0 and matlab 1
%     if ~isempty(location_idx) %) = location_idx;
%         IED_per_stimuli(spike,1) = location_idx; %Should we count all spikes?
%     elseif isempty(location_idx) %present in the intervals
%         IED_per_stimuli(spike,1) = nan;
%     end
% end %spike
