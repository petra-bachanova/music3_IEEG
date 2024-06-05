
function [parentfolder,path_to_toolboxes] = defineParentFolderAndPath_mus_DartFS()

% The function set the path to the required toolboxes and to the folder
% containing the analysis code. 
% ============================================================================
% *** Make sure the path below is correct before running the analysis code ***
% ============================================================================

clc;
parentfolder = '/Volumes/Ecog/music3_IEEG/Code/Matlab Scripts' ;    % path to the folder where the zip file was extracted
path_to_toolboxes =  '/Volumes/Ecog/music3_IEEG/Code/Matlab Scripts/Matlab_toolboxes';        % path to the folder containing the required MATLAB toolboxes (some of which were distributed with the analysis code)

fprintf('\n PATHS DEFINED. \n\n Parent folder: %s \n Toolboxed path: %s \n',parentfolder,path_to_toolboxes);

addpath(genpath(parentfolder)); %Adds paths to the current matlab session

end
