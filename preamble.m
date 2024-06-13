%% HOUSEKEEPING
clear
close all
clc
restoredefaultpath;

if any(computer== ["PCWIN64","MACI64"])  % Windows or Mac OS
    fileData.computerName = getenv('COMPUTERNAME');
    % Use defaults
    fileData.OutputPath = fullfile("Output files"); % Output directory
    fileData.WorkingPath = fullfile(pwd + "/"); % Working directory
elseif computer =="GLNXA64" % Linux OS
    fileData.computerName = getenv('COMPUTERNAME');
    %Use defaults
    fileData.OutputPath = fullfile("Output file"); % Output directory
    fileData.WorkingPath = fullfile(pwd + "\"); % Working directory
end

% Add required paths
addpath("atmospheric models\")
addpath("functions\")