function parseLogFile(filename)
    % Open the log file
    fid = fopen(filename, 'r');
    
    % Initialize variables to store paragraph text and table data
    paragraphText = {};
    tableData = {};
    
    % Read lines from the file
    tline = fgetl(fid);
    while ischar(tline)
        % Check if the line contains tabular data
        if contains(tline, ',')
            % Split the line into cells using comma delimiter
            rowData = strsplit(tline, ',');
            tableData = [tableData; rowData];
        else
            % Add the line to paragraph text
            paragraphText = [paragraphText; tline];
        end
        
        % Read the next line
        tline = fgetl(fid);
    end
    
    % Close the file
    fclose(fid);
% Display the results
disp('Paragraph Text:');
disp(paragraphText);
disp('Table Data:');
disp(tableData);

%% CLEAN DATA
%Extract EV markers
EV = extractBefore(paragraphText,'.wav -- Starttime');
EV = EV(~cellfun(@isempty,EV));
EV = extractAfter(EV, 'SoundEvent -- stim/newStim/');

%Add EV to table
idx_EV_behavior = find(contains(tableData(:,4), "Accuracy"));

%Sanity Check
if numel(idx_EV_behavior) ~= numel(EV)
    warning('Not all events have behavioral data associated with them')
    if numel(EV)>numel(idx_EV_behavior)
        warning(['The last event (' cell2mat(EV(5)) ') doesnt have behavioral data'])
    end    
end

%Append data to behavioral table
tableData{1,end+1} = 'EV';
tableData(end+1,1:end) = tableData(1,:); %padd


number_of_responses = unique(diff(idx_EV_behavior));
if numel(number_of_responses) ~=1
    warning('Number of responses per piece of music was different')
end
number_of_responses = number_of_responses-1; %Exclusive

for i = 1:numel(idx_EV_behavior)
    %Based on indexing
    %EV_names = repmat(EV(i),idx_EV_behavior(i+1) - idx_EV_behavior(i) -1 ,1);
    %tableData(idx_EV_behavior(i)+1: idx_EV_behavior(i+1) - 1,end) = EV_names;
    
    %Based on number of responses (~Should be constant!)
    EV_names = repmat(EV(i),number_of_responses ,1);
    tableData(idx_EV_behavior(i)+1: idx_EV_behavior(i)+number_of_responses,end) = EV_names;
end

%EVENTS
events_times = split(paragraphText,':');
events_times_table = cell2table(events_times, 'VariableNames',{'EV', 'start_time'});
%% CLEAN UP
%Remove repeated rows
tableData(idx_EV_behavior(2:end),:) = [];
%Remove padding
tableData(end,:) = [];

table = cell2table(tableData(2:end,:));
varNames = tableData(1,:);
table.Properties.VariableNames = varNames;

%Add music session
MusicSession = 1; %Fix me! Find by name
table.MusicSession = repmat(MusicSession, height(table),1); 
%Save
subjID = extractBetween(filename, '/mus', '/');
writetable(table, ['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/data/mus' cell2mat(subjID) '/ieeg/mus' cell2mat(subjID) '_SART_data.xlsx'])
writetable(events_times_table, ['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/data/mus' cell2mat(subjID) '/ieeg/mus' cell2mat(subjID) '_event_times.xlsx'])

end

%% FLAGS 04/26
%1.) Idenitfy song names and align them to the accuracy/reaction time
%2.) Append those names to the behavioral data
%3.) Extract toi from the file

