%% Calibration
% Script to run least squares regression for the load cell calibration. It
% will output 2 plots of raw data over time (Strain vs Time) and then the 
% plot of the linearized strain (Load vs. Strain)

clc
clear
close all

%% Import Data
data_file_path = "H:\Thesis\Load Cell\Multipoint Calibration\29cm\";
g0 = table2array(importfile2(data_file_path + "0g.txt", [1, Inf]));
g50 = table2array(importfile2(data_file_path + "50g.txt", [1, Inf]));
g100 = table2array(importfile2(data_file_path + "100g.txt", [1, Inf]));
g200 = table2array(importfile2(data_file_path + "200g.txt", [1, Inf]));
g300 = table2array(importfile2(data_file_path + "300g.txt", [1, Inf]));
g400 = table2array(importfile2(data_file_path + "400g.txt", [1, Inf]));
g500 = table2array(importfile2(data_file_path + "500g.txt", [1, Inf]));
g600 = table2array(importfile2(data_file_path + "600g.txt", [1, Inf]));
g700 = table2array(importfile2(data_file_path + "700g.txt", [1, Inf]));
g800 = table2array(importfile2(data_file_path + "800g.txt", [1, Inf]));
g900 = table2array(importfile2(data_file_path + "900g.txt", [1, Inf]));
g1000 = table2array(importfile2(data_file_path + "1000g.txt", [1, Inf]));



%% Plot Raw Data
Strain = vertcat(g0, g50, g100, g200, g300, g400, ...
    g500, g600, g700, g800, g900, g1000);
figure(1)
plot(Strain)
title("Raw Strain")
ylabel("Strain")

%% Force over Time
Average_zero = mean(g0);
Mass = vertcat(0 * ones(length(g0), 1), 50 * ones(length(g50), 1),...
    100 * ones(length(g100), 1), 200 * ones(length(g200), 1), ...
    300 * ones(length(g300), 1), 400 * ones(length(g400), 1), ...
    500 * ones(length(g500), 1), 600 * ones(length(g600), 1),...
    700 * ones(length(g700), 1), 800 * ones(length(g800), 1),...
    900 * ones(length(g900), 1), 1000 * ones(length(g1000), 1));
Force = (Mass) .*9.81./1000;
figure(2)
plot(Strain, Force)
title("Force vs. Strain")
xlabel("Strain")
ylabel("Force (N)")

%% Determine Slope of Line
% Perform least squares regression (linear regression, degree 1)
coefficients = polyfit(Strain, Force, 1);

% Extract the slope 'a' coefficient
a = coefficients(1);
b = coefficients(2);

%Add the function to the plot
functionText1 = sprintf('Least Squares Regression Function:');
functionText2 = sprintf('y = %.8fx + %.5f', a, b);
figure(2)
text(270000, 8.5, functionText1, 'FontSize', 12, 'HorizontalAlignment', 'left');
text(270000, 8, functionText2, 'FontSize', 12, 'HorizontalAlignment', 'left');
