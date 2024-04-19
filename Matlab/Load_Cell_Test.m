%% Calibration
% Script to run least squares regression for the load cell calibration. It
% will output 2 plots of raw data over time (Strain vs Time) and then the 
% plot of the linearized strain (Load vs. Strain)
% 
% clc
% clear
% close all
addpath("Import\");

%% Import Data
data_file_path = "H:\Thesis\Wave Tank Tests\Scaled ReTest\20240315-Perf2-2-0-1-Force.txt";
Strain = table2array(importfile2(data_file_path, [1, Inf]));



%% Plot Raw Data

figure(1)
plot(Strain(:))
title("Reference-Run #2-Test Wave #2-Raw Strain ")
xlabel("Bits")
ylabel("Strain")
x1 = [429 456 456 429];
x2 = [301 334 334 301];
y = [300000 300000 1200000 1200000];
% patch(x1,y,'red','FaceAlpha',.3);
% patch(x2,y,'green','FaceAlpha',.3);



%% Force over Time
Average_zero = mean(Strain(6:50));

%% Covert all data
Force = Strain(:) .* 0.000046040 - Average_zero*0.000046040;


%Plot Force Vs. Strain
% figure(2)
% plot(Force(50:end),Strain(50:end))
% title("Force vs. Strain")
% xlabel("Strain")
% ylabel("Force (N)")

%Plot Wave Force
figure()
plot(Force(50:end-1));
hold on
title("Wave Force")
xlabel("Readouts")
ylabel("Force (N)")

%%Max Wave Force
Force_Maxs(:) = maxk(abs(Force(72:100)),2,1);
Force_Max = mean(Force_Maxs);

[L0 L] = WaveSolver(AvgT(1), 0.4375, AvgH(1)*.01);
[Force1, Force_Maxs1, Force_Max1] = Load_Cell_Force(data_file_path,...
    6, 50,72, 100, L, AvgH(1)*.01);
%Plot Wave Force
figure()
plot(Force1);
hold on
title("Wave Force 1")
xlabel("Readouts")
ylabel("Force (N)")
