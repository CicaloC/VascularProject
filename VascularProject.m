clear all; close all; clc;
dataExp = load('Input_MarfanTAA_ATA.mat');
dataControl = load('Input_MarfanTAAcontrol_ATA.mat');
lambdaz_exp = dataExp.data_kl.lambdaz_exp;
Psys_exp = dataExp.data_kl.Psys_exp;
controlLambdaz_exp = dataExp.data_kl.lambdaz_exp;
controlPsys_exp = dataExp.data_kl.Psys_exp;