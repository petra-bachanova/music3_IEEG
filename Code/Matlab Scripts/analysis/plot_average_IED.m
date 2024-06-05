function [IEDs_morphology,x] = plot_average_IED(IEDs_Quon,data_eeg)

set_figure_colors
fs = data_eeg.hdr.Fs;
%% GET SPIKES ADJUSTED
IEDs_time_adjusted = [IEDs_Quon.spike_start_2048_hz IEDs_Quon.spike_end_2048_hz]; %Data comes in 2048 Hz
if ~(data_eeg.hdr.Fs == 2048)
    IEDs_time_adjusted = round((IEDs_time_adjusted/2048)*fs);
    %Explanation: spikes are recorded in full sampling rate (2048 Hz). If
    %the eeg file is not in 2048, adjust spikes to that sampling rate
    %(e.g., 2048 -> 512 Hz)
end


%% FIND SPIKES AND GET DATA IN ARRAY

IEDs_morphology = [];
IED_duration_start = fs*1; %how many seconds before/after start? Required to force vertcat
IED_duration_end = fs*3;
for spike_ID = 1:height(IEDs_Quon)
    current_chns = split(IEDs_Quon.channels(spike_ID),'_');
    current_chns = current_chns(~cellfun(@isempty,current_chns));
    for channel_ID = 1:numel(current_chns) %% REVISE!!!!!!
        %Ensure we take only IVs in the file
        if numel(data_eeg.time{1, 1})>(IEDs_time_adjusted(spike_ID,1)+IED_duration_end) %compare number of samples in file vs. end of file  
            IEDs_morphology = [IEDs_morphology; ... %ccat to previous morphology
            data_eeg.trial{1,1}(strcmp(data_eeg.label,current_chns(channel_ID)), ... %take data from the current channel
            IEDs_time_adjusted(spike_ID,1)-IED_duration_start:IEDs_time_adjusted(spike_ID,1)+IED_duration_end)]; %#ok<*AGROW>
        end
    end %channel_ID
end %spike_ID


%% PLOT AVERAGE IEDs
iv_highlight = [-0.5 1.00];

x = linspace(-IED_duration_start/fs,IED_duration_end/fs,width(IEDs_morphology));
y = mean(IEDs_morphology,1,"omitnan"); %data here is in eeg's fs.
h = figure('Name',['IED during experiment in one subject'],'position',[0 0 250 150],'color','w');
reduce_plot(x,y,"Color",COLOR.black); hold on;
reduce_plot(x(x >= iv_highlight(1) & x <= iv_highlight(2)), y(x >= iv_highlight(1) & x <= iv_highlight(2)),"Color",COLOR.red,'MarkerSize',1);
xlabel('Time from IED event (s)')
ylabel('voltage (μV)')
title(['Grand average in one subject. ' 'n = ' num2str(height(IEDs_morphology))])
ax = gca;
ax.FontSize = 16;
pause(2)

IEDs_singlesub = figure('Name',['IED during experiment in one subject'],'position',[0 0 250 150],'color','w'); hold on;
h1 = shadedErrorBar(x, mean(IEDs_morphology, "omitnan"), 1*std(IEDs_morphology,"omitmissing")./sqrt(size(IEDs_morphology,2)),{'color',COLOR.black,'lineWidth',1,'linesmoothing','on'},1);
alpha(h1.patch,0.25);
xlabel('Time from IED event (s)')
ylabel('voltage (μV)')
title(['Grand average in one subject. ' 'n = ' num2str(height(IEDs_morphology))])
xlim([-0.5 0.5])
ax = gca;
ax.FontSize = 16;
pause(2)
%% PLOT SAMPLE IEDS
iv_highlight = [-0.60 1.00];
h_examples = figure('Name',['IED sample traces'],'position',[0 0 500 500],'color','w'); hold on
num_examples = 12;
for example = 1:num_examples %modify as needed
    subplot(4,3,example);
    current_IED = IEDs_morphology(randi([0 height(IEDs_morphology)]),:); 
    reduce_plot(x,current_IED,"Color",COLOR.black,'MarkerSize',1); hold on;
    reduce_plot(x(x >= iv_highlight(1)& x <= iv_highlight(2)), current_IED(x >= iv_highlight(1)& x <= iv_highlight(2)),"Color",COLOR.red,'MarkerSize',1) 
    xlim([-1 2])
    axis off; axis tight;
end
