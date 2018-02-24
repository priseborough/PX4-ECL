clear all;
close all;

% add required paths
addpath('../Common');

% load compulsory data

dataIdentifier = 'APM_20170801';
baroFilePath = strcat('../TestData/',dataIdentifier,'/baro_data.mat');
load(baroFilePath);

imuFilePath = strcat('../TestData/',dataIdentifier,'/imu_data.mat');
load(imuFilePath);

magFilePath = strcat('../TestData/',dataIdentifier,'/mag_data.mat');
load(magFilePath);

gpsFilePath = strcat('../TestData/',dataIdentifier,'/gps_data.mat');
load(gpsFilePath);

% load optional data required for optical flow replay
rngFilePath = strcat('../TestData/',dataIdentifier,'/rng_data.mat');
flowFilePath = strcat('../TestData/',dataIdentifier,'/flow_data.mat');
if exist(rngFilePath,'file') && exist(flowFilePath,'file')
    load(rngFilePath);
    load(flowFilePath);
else 
    rng_data = [];
    flow_data = [];
end

% load optional data required for ZED camera replay
visoFilePath = strcat('../TestData/',dataIdentifier,'/viso_data.mat');
if exist(visoFilePath,'file')
    load(visoFilePath);
else
    viso_data = [];
end

% load default parameters
run('SetParameters.m');

%% run the filter replay 
output = RunFilter(param,imu_data,mag_data,baro_data,gps_data,rng_data,flow_data,viso_data);

%% generate and save output plots
runIdentifier = ' : APM data replay ';
folder = strcat('../OutputPlots/',dataIdentifier);
if ~exist(folder,'dir')
    mkdir(folder);
end
PlotData(output,folder,runIdentifier);

%% save output data
folder = strcat('../OutputData/',dataIdentifier);
fileName = strcat(folder,'/ekf_replay_output.mat');
if ~exist(folder,'dir')
    mkdir(folder);
end
save(fileName,'output');