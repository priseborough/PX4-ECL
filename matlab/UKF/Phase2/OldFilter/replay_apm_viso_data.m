clear all;
close all;

% add required paths
addpath('../Common');

% load compulsory data
load '../TestData/APM/baro_data.mat';
load '../TestData/APM/gps_data.mat';
load '../TestData/APM/imu_data.mat';
load '../TestData/APM/mag_data.mat';

rng_data = [];
flow_data = [];

% load data required for ZED camera replay
if exist('../TestData/APM/viso_data.mat','file')
    load '../TestData/APM/viso_data.mat';
else
    viso_data = [];
end

% load default parameters
run('SetParameters.m');

% ensure GPS is not used
param.control.gpsOnTime = 1e6;

% tuning of observation noise interpolation limits
%param.fusion.bodyVelErrorMin = 0.1; % Observation noise 1SD for the odometry sensor at the highest quality value (m/sec)
%param.fusion.bodyVelErrorMax = 0.9; % Observation noise 1SD for the odometry sensor at the lowest quality value (m/sec)

% run the filter replay 
output = RunFilter(param,imu_data,mag_data,baro_data,gps_data,rng_data,flow_data,viso_data);

% generate and save output plots
runIdentifier = ' : APM visual odometry replay ';
folder = strcat('../OutputPlots/APM_visual_odometry');
PlotData(output,folder,runIdentifier);

% save output data
folder = '../OutputData/APM_visual_odometry';
fileName = '../OutputData/APM_visual_odometry/ekf_replay_output.mat';
if ~exist(folder,'dir')
    mkdir(folder);
end
save(fileName,'output');