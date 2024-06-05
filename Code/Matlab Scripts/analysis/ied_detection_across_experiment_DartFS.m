%% IEEG-MUSIC: ANALYSIS OF IEDS DURING THE EXPERIMENT
%Author: CACR
%Date: 02/09/24
%[O]: Optional (do as required)
%[M]: Mandatory step

%BASICS
clc; clear;
addpath /Users/camilocastelblanco/Documents/GitHub/fieldtrip
ft_defaults
eegfile = 'Mus17_reclip_run1_BF_deidentified.edf'; %Δ
subjID = 'Mus17'; %Δ
MusicSession = 1; %Δ
plot_raw_spikes_flag = false;

cd('/Volumes/Ecog/music3_IEEG/Code/Matlab Scripts')
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath_mus_DartFS();

%LOAD BEHAVIORAL DATA
trial_data = readtable("/Volumes/Ecog/music3_IEEG/Data/music3_trialdata_mc.xlsx", Sheet="Main");
SART_data = readtable("/Volumes/Ecog/music3_IEEG/Data/music3_SARTdata.xlsx", Sheet="Main");

% LOAD  iEEG DATA
eegfile = ['/Volumes/Ecog/music3_IEEG/Data/EDF_files/', subjID,'_reclip_run1_BF_deidentified.edf'];
hdr = ft_read_header(eegfile);
fs = hdr.Fs;

%% PREPROCESS LOG DATA [O/M]
parseLogFile(['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/data/' subjID '/ieeg/17_AK.log']);

%% COMPUTE ACCURACY/REACTION TIME 
process_SART_data(['/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/data/' subjID '/ieeg/' subjID '_SART_data.xlsx']);

%% LOAD IEDs Quon 2022 [M]
cd('/Volumes/Ecog/music3_IEEG/Code/Matlab Scripts/preprocessing') %Unsure how to add path in cluster
cfg_open_IEDs = [];
cfg_open_IEDs.subjID = subjID; %Defined above
cfg_open_IEDs.downsample_fs = 200; %sampling rate used in the CNN (Check ipynb)
cfg_open_IEDs.recording = 'run1'; %if only one recording, leave empty. otherwise, specify r1/r2
cfg_out_IEDs = import_IEDs_Quon2022_csv_mus_DartFS(cfg_open_IEDs);

%% LOAD DATA ON EVENT MARKERS [M]
cfg_IEDs_events = struct;
cfg_IEDs_events.IEDs_Quon = cfg_out_IEDs; 
cfg_IEDs_events.trial_data = trial_data;
cfg_IEDs_events.SART_data = SART_data;
cfg_IEDs_events.subjID = subjID;
cfg_IEDs_events.MusicSession = MusicSession;
cfg_IEDs_events.hdr = hdr;
cfg_IEDs_events.rate_of_interest = 30;
cfg_IEDs_events.bin_time = 15;
cfg_out_IEDs_with_events = ied_detection_on_event_markers_DartFS(cfg_IEDs_events, subjID);

%% SAVE
cfg_out_IEDs_with_events.channels_frequency = cfg_out_IEDs.channels_frequency;  
cfg_out_IEDs_with_events.chn_dictionary = cfg_out_IEDs.chn_dictionary;  
save(strcat('/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/scripts/Matlab Scripts/preprocessing/results/preprocessed_IED_rates/IEDs_',subjID, '_MusicSession_', num2str(MusicSession)), "cfg_out_IEDs_with_events");

%% PREPROCESS iEEG DATA TO PLOT [O]
cfg = [];
cfg.dataset = eegfile;
eegchan          = strcat('-', ft_channelselection({'eeg'}, hdr.label));
cfg.channel      = ft_channelselection({'all', '-*DC*', '-PR', '-Pleth', '-TRIG', '-OSAT', '-C2*', '-C1*',eegchan{:}}, hdr.label); %#ok<CCAT>

if plot_raw_spikes_flag
    data_eeg    = ft_preprocessing(cfg); %#ok<*UNRCH>
    cfg = [];
    cfg.bsfilter     = 'yes'; %https://www.allaboutcircuits.com/uploads/articles/6.13_Band-Pass_and_Band-Reject_Active_Filters1_.jpg
    cfg.bsfiltord    = 3; 
    cfg.bsfilttype   = 'but'; %default
    cfg.bsfiltdir = 'twopass'; %filter direction
    cfg.bsfiltwintype = 'hamming';
    cfg.bsfreq       = [59 61; 119 121; 179 181];
    cfg.lpfilter     = 'yes';
    cfg.lpfreq       = 250;
    data_eeg_filtered = ft_preprocessing(cfg, data_eeg);
    %data_eeg_filtered = ft_preprocessing(cfg, data_eeg_filtered);
    
    if fs > 512
        cfg = [];
        cfg.resamplefs = 512;
        cfg.method = 'resample';
        %data_eeg_filtered = ft_resampledata(cfg, data_eeg_filtered);
        data_eeg_filtered = ft_resampledata(cfg,data_eeg);
    end
end

% PLOT AVERAGE IEDs
IEDs_Quon = cfg_out_IEDs_with_events.IED_table;
[average_IED,x] = plot_average_IED(IEDs_Quon,data_eeg_filtered); %THE IEDs_Quon TABLE IS IN 2048 FS

%% PSD ANALYSIS
[S,freq] = spectopo(average_IED,frames,fs,'freqrange',[2 floor(fs/2)],'percent',50,'plot','off');

figure(); hold on;
set(findobj(gca,'type','line'),'Color',[0.5 0.5 0.5])
frames = 0; %	frames per epoch {0 -> whole data length}
PSD_lines_fig = plot(freq,mean(S,1,"omitnan"),'k'); %mean(S,1,"omitnan")
ylim_S = min(S, [],'all');
xlim([0 200]); ylim([ylim_S 50]);
xlabel('Frequencies (Hz)'); ylabel('dB');
title('Power Spectral Density')
ax = gca; 
ax.FontSize = 16;

save('power_spectrum', 'freq', 'S');


%% PREPROCESS DATA

%FIELDTRIP STRUCTURES
data.label = {'avg'};     % cell-array containing strings, Nchan*1
data.fsample = fs;   % sampling frequency in Hz, single number
data.trial =  {mean(average_IED,1,"omitnan")};   % cell-array containing a data matrix for each
                % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
data.time =  {x}; % cell-array containing a time axis for each
                % trial (1*Ntrial), each time axis is a 1*Nsamples vector


%% TIME FREQUENCY ANALYSIS
time_res = 0.002;
cfg = []; 
cfg.method = 'wavelet'; %mtmfft
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.foi = [2:2:300];
cfg.toi = data.time{1}(1):time_res:data.time{1}(end);
cfg.pad   = 'nextpow2';
cfg.polyremoval = 0; %equivalent to demeaning
cfg.width  = 7; %number of cycles, of the wavelet
cfg.gwidth = 3; %determines the length of the used wavelets in standard deviations
freq = ft_freqanalysis(cfg, data);

%% BASELINE NORMALIZATION           
cfg = [];
cfg.baseline = [-.7 -.5]; %-0.8 -0.6
cfg.baselinetype = 'absolute'; %see line 196 (https://github.com/fieldtrip/fieldtrip/blob/master/ft_freqbaseline.m)
freq_blc = ft_freqbaseline(cfg, freq); 

%% PLOT SPECTROGRAM

mean_specto = figure('Name',['Mean IED'],'position',[0 0 1500 500]);
imagesc(freq_blc.time, freq_blc. freq, squeeze(freq_blc.powspctrm), [-0.5 20]); %last parameter delimiters z axis
axis xy % flip vertically
xlabel("Time from IED spike start(s)"); ylabel("Frequency (Hz)");
title('Mean IED');
ylim([6 200]);
xlim([-0.5 1]); 
ax = gca;
ax.FontSize = 16;
colorbar
set(gca,'fontname','SansSerif') 
set(gcf,'color','w') %Remove ugly white background

%%
x = linspace(1,10,100); a = exp(1); %Euler
y_viewscrambled = linspace(0,length(x)*0.5,length(x));
y_viewnonfamiliar = linspace(0,length(x)*0.75,length(x));
%SIGMOID NEURON: 
%y = 1/(1+e^(wx +b))
w = -2.25; %"sharpness" of sigmoid
b = 17; %x shifting
y_familiar = length(x)./(1+a.^(w.*x+b));
%https://prvnk10.medium.com/sigmoid-neuron-ad0ec6f9a3e2
%https://www.desmos.com/calculator

figure(); hold on;
plot(x,y_viewscrambled, 'Linewidth', 5)
plot(x,y_viewnonfamiliar, 'Linewidth', 5)
plot(x,y_familiar, 'Linewidth', 5)
mycolors = [COLOR.blue; COLOR.green; COLOR.reddark];
set(gca,'ColorOrder',mycolors);

xlabel('Time');  xticklabels([""]); yticks([])
ylabel('Neocortical Activity'); yticklabels([""]);xticks([])
axis tight;

legend({'Scrambled','Non-familiar', 'Familiar'});
legend boxoff
set(legend,'fontsize',16,'Location','NorthWest','LineWidth',1)

memory_recall = 6;
valence_arrow = annotation('textarrow'); valence_arrow.Parent = gca;
valence_arrow.X = [memory_recall+1,memory_recall]; % set the x-property
valence_arrow.Y = [2,2];
valence_arrow.String = 'Memory Recall';
valence_arrow.FontSize = 14;
ax = gca;
ax.FontSize = 16;
set(gcf,'color','w'); %Set background to white

xline(memory_recall, '--', 'HandleVisibility','off')