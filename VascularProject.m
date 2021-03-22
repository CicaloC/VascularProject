clear all; close all; clc;
%%
global inputDataExp inputDataControl

inputDataExp = load('Input_MarfanTAA_ATA.mat');
inputDataControl = load('Input_MarfanTAAcontrol_ATA.mat');


% Pressure & Axial Stretch @ loaded configurations of interest
lambdaz = inputDataExp.data_kl.lambdaz_exp; %stored in input data

Pi = inputDataExp.data_kl.Psys_exp; %in-vivo pressure

%% Solve the equilibrium equation
% Equilibrium equation for a pressurized and axially stretched vessel
% Unknown: outer radius, ro
% ro will be estimated as the roots of the equilibrium equation

H=cell(1,1);
H{1}=@equilibrium_r_or_loaded; %handle to equilibrium equation 1
H{2}=@equilibrium_z_fz_loaded; %handle to equilibrium equation 2

% Parameter to be estimated in loaded configuration is the outer radius (ro) 
% x0 is the initial guess for ro based on the input experimental data
x0 = [inputDataExp.data_ktf.or_exp];

       
inputDataExp.data_kl.lambdaz = lambdaz; % axial stretch
inputDataExp.data_kl.Pi = Pi; % liminal pressure, in kPa
        
% Newton Raphson Scheme, Inputs: function handles and evaluation points
x = Newton_Raphson(H{2},lambdaz); % call the Newton Raphson function (Newton_Raphson_tutorial.m)
or_est = x(1,1); % calculated outer radius
        
