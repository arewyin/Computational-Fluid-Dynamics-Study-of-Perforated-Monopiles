%% Load Cell Readout
% Script reads out a singular strain set and then converts to force using
% calibration slope

clc
clear
close all

%% Import Data
Strain_Ref_8 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run1-7_4hz-REF-Force.txt", [1, Inf]));
Strain_Perf1_8_0 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run2-7_4hz-Perf1_0deg-Force.txt", [1, Inf]));
Strain_Perf1_8_15 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run1-7_4hz-Perf1_15deg-Force.txt", [1, Inf]));
Strain_Perf1_8_30 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run1-7_4hz-Perf1_30deg-Force.txt", [1, Inf]));
Strain_Perf1_8_45 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run2-7_4hz-Perf1_45deg-Force.txt", [1, Inf]));

Strain_Perf2_8_0 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run1-7_4hz-Perf2_0deg-Force.txt", [1, Inf]));
Strain_Perf2_8_45 = table2array(importfile2("H:\Thesis\Wave Tank Tests\Scaled Tests\20231117-Run1-7_4hz-Perf2_45deg-Force.txt", [1, Inf]));

%% Plot Raw Data
figure(1)
plot(Strain_Perf1_8_45)
title("Raw Strain vs. Time")
xlabel("Readouts")
ylabel("Strain")

%% Zero the data
Zeroed_Strain_Ref_8 = mean(Strain_Ref_8(786:end-1));
Zeroed_Strain_Perf1_8_0 = mean(Strain_Perf1_8_0(6:43));
Zeroed_Strain_Perf1_8_15 = mean(Strain_Perf1_8_15(1:73));
Zeroed_Strain_Perf1_8_30 = mean(Strain_Perf1_8_30(1:76));
Zeroed_Strain_Perf1_8_45 = mean(Strain_Perf1_8_45(1:76));

Zeroed_Strain_Perf2_8_0 = mean(Strain_Perf2_8_0(1:20));
Zeroed_Strain_Perf2_8_45 = mean(Strain_Perf2_8_45(99:end));

%% Convert strain to force
Force_Ref_8 = Strain_Ref_8(786:end-1) .* 0.0008 - Zeroed_Strain_Ref_8*0.0008;
Force_Perf1_8_0 = Strain_Perf1_8_0(6:end-1) .* 0.0008 - Zeroed_Strain_Perf1_8_0*0.0008;
Force_Perf1_8_15 = Strain_Perf1_8_15(1:end-1) .* 0.0008 - Zeroed_Strain_Perf1_8_15*0.0008;
Force_Perf1_8_30 = Strain_Perf1_8_30(1:end-1) .* 0.0008 - Zeroed_Strain_Perf1_8_30*0.0008;
Force_Perf1_8_45 = Strain_Perf1_8_45(1:end-1) .* 0.0008 - Zeroed_Strain_Perf1_8_45*0.0008;

Force_Perf2_8_0 = Strain_Perf2_8_0(6:end-1) .* 0.0008 - Zeroed_Strain_Perf2_8_0*0.0008;
Force_Perf2_8_45 = Strain_Perf2_8_45(99:end) .* 0.0008 - Zeroed_Strain_Perf2_8_45*0.0008;

%% Replot data
figure(2)
plot(Force_Ref_8)
hold on
% plot(Force_Perf1_8_0)
% plot(Force_Perf1_8_15)
% plot(Force_Perf1_8_30)
% plot(Force_Perf1_8_45)
% plot(Force_Perf2_8_0)
plot(Force_Perf2_8_45)
title("Force vs. Time")
xlabel("Readouts")
ylabel("Force")
% legend('Reference','Perforation 1 - Deg 0','Perforation 1 - Deg 15',...
%     'Perforation 1 - Deg 30','Perforation 1 - Deg 45','Perforation 2 - 0 Deg'...
%     ,'Perforation 2 - 45 Deg')

%% Calculate average force max
Force_Max_Ref_8 = mean(maxk(abs(Force_Ref_8),20,1));
Force_Max_Perf1_8_0 = mean(maxk(abs(Force_Perf1_8_0),20,1));
Force_Max_Perf1_8_15 = mean(maxk(abs(Force_Perf1_8_15),20,1));
Force_Max_Perf1_8_30 = mean(maxk(abs(Force_Perf1_8_30),20,1));
Force_Max_Perf1_8_45 = mean(maxk(abs(Force_Perf1_8_45),20,1));
Force_Max_Perf2_8_0 = mean(maxk(abs(Force_Perf2_8_0),20,1));
Force_Max_Perf2_8_45 = mean(maxk(abs(Force_Perf2_8_45),20,1));