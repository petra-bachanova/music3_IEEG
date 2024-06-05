%% CONVERTING FINAL SPIKES FROM QUON 2022 INTO MATLAB
function cfg_out_IEDs = import_IEDs_Quon2022_csv_mus_DartFS(cfg_open_IEDs)

%% LOAD DATA
subjID = cfg_open_IEDs.subjID;
downsample_fs = cfg_open_IEDs.downsample_fs;
recording = cfg_open_IEDs.recording;
spike_data = readtable(['/Volumes/Ecog/music3_IEEG/Data/Detected_spikes/',subjID ,'_', recording,'_finalspikes_240603.xlsx']);
%channelsData = readtable(['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/scripts/Matlab Scripts/IEDs_Quon2022/',subjID, '/',subjID,recording,'_finalspikes_channels.xlsx']);

%% GET DATA
%Spike start in dsampled data
spike_start_Quon = spike_data.spike_start;
%Channels
chn_dictionary = dictionary(spike_data.channels_idx, spike_data.channels_names); %create dictionary
strToRemove = ['"', "'", "[", "]", ','];
spike_data.channels = erase(spike_data.channels_names,strToRemove);
channel_numbers = spike_data.channels_idx;
channels = spike_data.channels_names;

channels_wordcloud = countWordsInArray(spike_data.channels);
channels_wordcloud(1,:) = [];
%%
% for allSpikes = 1:height(spike_data)
%     current_channels = cell2mat(spike_data.channels(allSpikes)); %get channels in current spike
%     current_channels = split(current_channels,' '); %separate based on spaces
%     %current_channels = str2double(current_channels);
%     chn_str = '';
%     for chn = 1:numel(current_channels)
%         chn_str = strcat(chn_str,cell2mat(chn_dictionary(current_channels(chn,1))),'_');
%         channels_wordcloud = [channels_wordcloud; chn_dictionary(current_channels(chn,1))]; %#ok<AGROW>
%     end %chn
%     channels(allSpikes,1) = {chn_str};
% end %allSpikes 
% clear strToRemove allSpikes channelsData chn_str current_channels strToRemove

%% CONVERT FROM DS DATA TO FULLY SAMPLED DATA
IED_duration = 0.05*2048; %remove samples after start

spike_start_2048_hz = convert_ivs_to_other_fs(spike_data.spike_start, downsample_fs, 2048);
spike_end_2048_hz = round(spike_start_2048_hz+IED_duration);
num_spikes = spike_data.channels_count; %gets num of IEDs from the number of channels detected
spike_start_Quon_CNN = spike_data.spike_start;
spike_start_seconds = spike_start_2048_hz/2048;
IEDs = table(spike_start_Quon_CNN,spike_start_2048_hz,spike_end_2048_hz,spike_start_seconds,channels,channel_numbers,num_spikes);


%% ANALYSIS OF THE LOCATIONS WHERE IEDs ARE FOUND
% uniqueChn = unique(channels_wordcloud);
% counter = [];
% for chn = uniqueChn' %needs to be transposed to iterate
%     counter = [counter ; sum(strcmp(channels_wordcloud, chn))]; %#ok<AGROW>
% end
% channels_frequency = table(uniqueChn, counter);
figure()
wordcloud(channels_wordcloud,'Word','Count');
set(gca,'fontname','SansSerif') 
set(gcf,'color','w') 
title('IEDs: Detection (Quon 2022). Electrode locations)')
saveas(gcf, ['/Volumes/Ecog/music3_IEEG/Results/spike_wordcloud/', subjID, '_worcloud.jpg'])
pause(0.3)
close
%% SAVE OUTPUT
cfg_out_IEDs = struct;
cfg_out_IEDs.IEDs = IEDs;
cfg_out_IEDs.channel_wordcloud = channels_wordcloud;
% cfg_out_IEDs.channels_frequency = channels_frequency;
% cfg_out_IEDs.chn_dictionary = chn_dictionary;
cfg_out_IEDs.subject = subjID;

end

function wordCounts = countWordsInArray(arr)
    % Initialize a container for word counts
    wordMap = containers.Map();
    
    % Loop through each row of the array
    for i = 1:size(arr, 1)
        % Split the row into words using space as delimiter
        words = strsplit(arr{i});
        
        % Loop through each word and update its count in the map
        for j = 1:length(words)
            word = words{j};
            if isKey(wordMap, word)
                wordMap(word) = wordMap(word) + 1;
            else
                wordMap(word) = 1;
            end
        end
    end
    
    % Convert the map to a table for better readability
    keys = wordMap.keys();
    values = wordMap.values();
    wordCounts = table(keys', [values{:}]', 'VariableNames', {'Word', 'Count'});
end

%% CODE VAULT
% %How many spikes are there?
% numChannels = spikesData.numChannels;
% totalNumSpikes = sum(numChannels); clear numChannels
% 
% %Predefine matrix
% dataSorted = zeros(totalNumSpikes, 4); 
% % 1st row = channel; 2nd = start, duration is unkown
% 
% %How many events (i.e., spikes) we have ?
% numOfEvents =  height(spikesData);
% 
% %remove '' from table (working on making automatic... for now done on excel
% strToRemove = ['"', "'", "[", "]", ','];
% spikesData.channels = erase(spikesData.channels,strToRemove);
% 
% %Don't forget remove columns that aren't clean
% %% Sorting out data
% 
% save_flag = 0; %set to 1 to save output
% dataSortedCounter = 1; %Note loop does not account for table size as it iterates differently due to data nature
% for i = 1:numOfEvents
%     str = cell2mat(spikesData.channels(i)); %get channels in current spike
%     channels = split(str,' '); %separate based on spaces
%     channels = str2double(channels);
%     for ii = 1:length(channels)
%         dataSorted(dataSortedCounter, 1) = channels(ii); %channel
%         dataSorted(dataSortedCounter, 2) = spikesData.spikeStart(i); %spike start
%         dataSorted(dataSortedCounter, 3) = char(i); %idx
%         dataSorted(dataSortedCounter, 4) = length(channels); %number of spikes per
%         %detection
%         dataSortedCounter = dataSortedCounter+1;
%     end
% end
% channel = dataSorted(:,1);
% start = dataSorted(:,2);
% idx = dataSorted(:,3);
% numOfSpikesperDetection  = dataSorted(:,4);
% dataSortedTable = table(channel, start, idx, numOfSpikesperDetection);
% 
% 
% if save_flag == 1
%     save('sub-01_IEDsDouble.mat', "dataSorted");
%     save('sub-01_EDsTable.mat', "dataSortedTable");
% end
% %% Distrbution of spikes in time
% 
% tableToPlot = zeros(max(dataSorted(:,2)),1); %create full timeline
% Times = seconds(linspace(1,max(dataSorted(:,2)), max(dataSorted(:,2)))');
% for aa = 1:numOfEvents %go through all spike detections
%     idx_time = spikesData.spikeStart(aa); %find its time
%     tableToPlot(idx_time,1) = spikesData.numChannels(aa); %add freques to table
% end
% tableToPlot = array2table(tableToPlot, 'VariableNames',{'Frequency'});
% dataSortedTimeTable = table2timetable(tableToPlot, 'RowTimes', Times);
% s = stackedplot(dataSortedTimeTable,{'Frequency'}, 'Title', [{'sub-01: IEDs'}],'DisplayLabels',{'Frequency'});
% set(gca,'fontsize', 14)
% 
% %% Distribution of spikes in time (initial)
% 
% %discontinuos data
% Time = seconds(dataSorted(:,2)); %create time vector
% dataSortedTimeTable = table2timetable(dataSortedTable, 'RowTimes',Time);
% s = stackedplot(dataSortedTimeTable,{'numOfSpikesperDetection'}, 'Title', 'sub-01_IEDs','DisplayLabels',{'Number of Spikes'});
% s.LineProperties(1).PlotType = 'stairs';
% s.XLabel = {'Time(s)'}
%%
% writematrix('data', 'behavioralDataPreprocessed.xlsx')
% writetable(data,'behavioralDataPreprocessed.xlsx')