%TITLE: mus_28 preprocessing
%AUTHOR: Camilo Castelblanco Riveros % Brian Fidali
%DATE: 12/07/23


clc; clear;
addpath /Users/camilocastelblanco/Documents/GitHub/fieldtrip
ft_defaults
eegfile = 'mus28.EDF';
set_figure_colors
subjID = 'mus_28';
cd( '/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/Scripts/Matlab Scripts')
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath_mus();

%% LOAD AND PROCESS BEHAVIORAL DATA 
cd('/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/mus28_test_20231116')
behevaioral_data = readtable("subject28_13.55.log", 'FileType','text');
test = readmatrix("subject28_13.55.log", 'Delimiter',{',','='},'Range','A1');
test2 =  num2cell(fscanf(fopen("subject28_13.55.log"), '%f')');
%% DEFINE TRIALS
cd( '/Users/camilocastelblanco/Dropbox (Dartmouth College)/ieeg-music/mus28_test_20231116/ieeg')
cfg              = [];
cfg.dataset      = eegfile;
cfg.trialfun = 'mus28_trialDefinition';
cfg.trialdef.pre  = 1; %before flanking we shall take 1 seconds
cfg.trialdef.post = 1;
cfg = ft_definetrial(cfg);

