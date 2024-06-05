%% IEEG-MUSIC: ANALYSIS OF ERROS
%Author: CACR
%Date: 02/09/24

clc; clear;
addpath /Users/camilocastelblanco/Documents/GitHub/fieldtrip
ft_defaults
eegfile = 'Subject18_Session 1_deidentified.edf';
set_figure_colors
subjID = 'mus18';
MusicSession = 1;
cd( '/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/Scripts/Matlab Scripts')
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath_mus();

%%
cd('/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/scripts/Python Scripts/MUS18_bcf')
eeg_data = readtable('mus18_r1_SEEG+scalp_eegdata.csv');
eeg_data_array = table2array(eeg_data);
eeg_data_nan = find(isnan(eeg_data_array));