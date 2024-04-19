%% Scaled Monopile Data Analysis
% For the data analysis of forces on scaled monopiles with and without
% perforations. This script will output 3 graphs for each of the
% perforation geometries relating the force on the perforation/force on the
% reference vs. KC number. Each graph contains the trends for 4 wave
% approach angles (0, 15, 30, and 45) at 3 wave test states (6, 8, 13) and
% the applicable trend from either Anderson, et al. (2020) and Ploeg (2021)

clc
clear
close all
addpath("Import\");

%%
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontSize', 12) 
set(0,'defaultAxesFontName','Arial')

%% Variable Setup
KC = [1.16, 5.15, 2.12];
datafilepath = "H:\Thesis\Scaled Simulations\Scaled Data Analysis\";

%% Data Import
% Reference (3 total)
Ref_6 = table2array(import01(datafilepath + "Reference_6.txt", 9, inf));
Ref_8 =  table2array(importPressAndForce(datafilepath + "Reference_8.txt", 9, inf));
Ref_13 =  table2array(import01(datafilepath + "Reference_13.txt", 9, inf));

% Perforation 1 (12 total)
Perf_1_6_0 =  table2array(import01(datafilepath + "Perforation_1_6_0.txt", 9, inf));
Perf_1_6_15 = table2array(import01(datafilepath + "Perforation_1_6_15.txt", 9, inf));
Perf_1_6_30 = table2array(import01(datafilepath + "Perforation_1_6_30.txt", 9, inf));
Perf_1_6_45 = table2array(import01(datafilepath + "Perforation_1_6_45.txt", 9, inf));
Perf_1_8_0 = table2array(import01(datafilepath + "Perforation_1_8_0_New.txt", 9, inf));
Perf_1_8_15 = table2array(import01(datafilepath + "Perforation_1_8_15_New.txt", 9, inf));
Perf_1_8_30 = table2array(import01(datafilepath + "Perforation_1_8_30_New.txt", 9, inf));
Perf_1_8_45 = table2array(import01(datafilepath + "Perforation_1_8_45_New.txt", 9, inf));
Perf_1_13_0 = table2array(import01(datafilepath + "Perforation_1_13_0.txt", 9, inf));
Perf_1_13_15 = table2array(import01(datafilepath + "Perforation_1_13_15.txt", 9, inf));
Perf_1_13_30 = table2array(import01(datafilepath + "Perforation_1_13_30.txt", 9, inf));
Perf_1_13_45 = table2array(import01(datafilepath + "Perforation_1_13_45.txt", 9, inf));

% Perforation 2 (12 total)
Perf_2_6_0 = table2array(import01(datafilepath + "Perforation_2_6_0.txt", 9, inf));
Perf_2_6_15 = table2array(import01(datafilepath + "Perforation_2_6_15.txt", 9, inf));
Perf_2_6_30 = table2array(import01(datafilepath + "Perforation_2_6_30.txt", 9, inf));
Perf_2_6_45 = table2array(import01(datafilepath + "Perforation_2_6_45.txt", 9, inf));
Perf_2_8_0 = table2array(import01(datafilepath + "Perforation_2_8_0_new.txt", 9, inf));
Perf_2_8_15 = table2array(import01(datafilepath + "Perforation_2_8_15_New.txt", 9, inf));
Perf_2_8_30 = table2array(import01(datafilepath + "Perforation_2_8_30_New.txt", 9, inf));
Perf_2_8_45 = table2array(import01(datafilepath + "Perforation_2_8_45_New.txt", 9, inf));
Perf_2_13_0 = table2array(import01(datafilepath + "Perforation_2_13_0.txt", 9, inf));
Perf_2_13_15 = table2array(import01(datafilepath + "Perforation_2_13_15.txt", 9, inf));
Perf_2_13_30 = table2array(import01(datafilepath + "Perforation_2_13_30.txt", 9, inf));
Perf_2_13_45 = table2array(import01(datafilepath + "Perforation_2_13_45.txt", 9, inf));

% Perforation 3 (12 total).
Perf_3_6_0 = table2array(import01(datafilepath + "Perforation_3_6_0.txt", 9, inf));
Perf_3_6_15 = table2array(import01(datafilepath + "Perforation_3_6_15.txt", 9, inf));
Perf_3_6_30 = table2array(import01(datafilepath + "Perforation_3_6_30.txt", 9, inf));
Perf_3_6_45 = table2array(import01(datafilepath + "Perforation_3_6_45.txt", 9, inf));
Perf_3_8_0 = table2array(import01(datafilepath + "Perforation_3_8_0_New.txt", 9, inf));
Perf_3_8_15 = table2array(import01(datafilepath + "Perforation_3_8_15_New.txt", 9, inf));
Perf_3_8_30 = table2array(import01(datafilepath + "Perforation_3_8_30_New.txt", 9, inf));
Perf_3_8_45 = table2array(import01(datafilepath + "Perforation_3_8_45_New.txt", 9, inf));
Perf_3_13_0 = table2array(import01(datafilepath + "Perforation_3_13_0.txt", 9, inf));
Perf_3_13_15 = table2array(import01(datafilepath + "Perforation_3_13_15.txt", 9, inf));
Perf_3_13_30 = table2array(import01(datafilepath + "Perforation_3_13_30.txt", 9, inf));
Perf_3_13_45 = table2array(import01(datafilepath + "Perforation_3_13_45.txt", 9, inf));


%% Average of Peak Total Force
% Reference
Ref_6_Peak = mean(maxk(Ref_6(4:150,3),10,1));
Ref_8_Peak = mean(maxk(Ref_8(3:750,4),6,1));
Ref_13_Peak = mean(maxk(Ref_13(3:150,3),6,1));

% Perforation 1
Perf_1_6_0_Peak = mean(maxk(Perf_1_6_0(3:200,3),6,1));
Perf_1_6_15_Peak = mean(maxk(Perf_1_6_15(3:200,3),6,1));
Perf_1_6_30_Peak = mean(maxk(Perf_1_6_30(3:200,3),6,1));
Perf_1_6_45_Peak = mean(maxk(Perf_1_6_45(3:200,3),6,1));
Perf_1_8_0_Peak = mean(maxk(Perf_1_8_0(3:530,3),6,1));
Perf_1_8_15_Peak = mean(maxk(Perf_1_8_15(3:end,3),6,1));
Perf_1_8_30_Peak = mean(maxk(Perf_1_8_30(3:end,3),6,1));
Perf_1_8_45_Peak = mean(maxk(Perf_1_8_45(3:530,3),6,1));
Perf_1_13_0_Peak = mean(maxk(Perf_1_13_0(3:150,3),6,1));
Perf_1_13_15_Peak = mean(maxk(Perf_1_13_15(3:150,3),6,1));
Perf_1_13_30_Peak = mean(maxk(Perf_1_13_30(3:150,3),6,1));
Perf_1_13_45_Peak = mean(maxk(Perf_1_13_45(3:150,3),6,1));

% Perfotation 2
Perf_2_6_0_Peak = mean(maxk(Perf_2_6_0(3:200,3),6,1));
Perf_2_6_15_Peak = mean(maxk(Perf_2_6_15(3:200,3),6,1));
Perf_2_6_30_Peak = mean(maxk(Perf_2_6_30(3:200,3),6,1));
Perf_2_6_45_Peak = mean(maxk(Perf_2_6_45(3:200,3),6,1));
Perf_2_8_0_Peak = mean(maxk(Perf_2_8_0(3:end,3),6,1));
Perf_2_8_15_Peak = mean(maxk(Perf_2_8_15(3:end,3),6,1));
Perf_2_8_30_Peak = mean(maxk(Perf_2_8_30(3:end,3),6,1));
Perf_2_8_45_Peak = mean(maxk(Perf_2_8_45(3:530,3),6,1));
Perf_2_13_0_Peak = mean(maxk(Perf_2_13_0(3:150,3),6,1));
Perf_2_13_15_Peak = mean(maxk(Perf_2_13_15(3:150,3),6,1));
Perf_2_13_30_Peak = mean(maxk(Perf_2_13_30(3:150,3),6,1));
Perf_2_13_45_Peak = mean(maxk(Perf_2_13_45(3:150,3),6,1));

% Perforation 3
Perf_3_6_0_Peak = mean(maxk(Perf_3_6_0(3:200,3),6,1));
Perf_3_6_15_Peak = mean(maxk(Perf_3_6_15(3:200,3),6,1));
Perf_3_6_30_Peak = mean(maxk(Perf_3_6_30(3:200,3),6,1));
Perf_3_6_45_Peak = mean(maxk(Perf_3_6_45(3:200,3),6,1));
Perf_3_8_0_Peak = mean(maxk(Perf_3_8_0(3:end,3),6,1));
Perf_3_8_15_Peak = mean(maxk(Perf_3_8_15(3:533,3),6,1));
Perf_3_8_30_Peak = mean(maxk(Perf_3_8_30(3:533,3),6,1));
Perf_3_8_45_Peak = mean(maxk(Perf_3_8_45(3:533,3),6,1));
Perf_3_13_0_Peak = mean(maxk(Perf_3_13_0(3:150,3),6,1));
Perf_3_13_15_Peak = mean(maxk(Perf_3_13_15(3:150,3),6,1));
Perf_3_13_30_Peak = mean(maxk(Perf_3_13_30(3:150,3),6,1));
Perf_3_13_45_Peak = mean(maxk(Perf_3_13_45(3:150,3),6,1));

%% Relation between Perforation Force and Total Force
% Perforation 1
Rel_F_1_0 = [Perf_1_6_0_Peak/Ref_6_Peak, Perf_1_8_0_Peak/Ref_8_Peak, ...
    Perf_1_13_0_Peak/Ref_13_Peak];
Rel_F_1_15 = [Perf_1_6_15_Peak/Ref_6_Peak, Perf_1_8_15_Peak/Ref_8_Peak, ...
    Perf_1_13_15_Peak/Ref_13_Peak];
Rel_F_1_30 = [Perf_1_6_30_Peak/Ref_6_Peak, Perf_1_8_30_Peak/Ref_8_Peak, ...
    Perf_1_13_30_Peak/Ref_13_Peak];
Rel_F_1_45 = [Perf_1_6_45_Peak/Ref_6_Peak, Perf_1_8_45_Peak/Ref_8_Peak, ...
    Perf_1_13_45_Peak/Ref_13_Peak];

% Perforation 2
Rel_F_2_0 = [Perf_2_6_0_Peak/Ref_6_Peak, Perf_2_8_0_Peak/Ref_8_Peak, ...
    Perf_2_13_0_Peak/Ref_13_Peak];
Rel_F_2_15 = [Perf_2_6_15_Peak/Ref_6_Peak, Perf_2_8_15_Peak/Ref_8_Peak, ...
    Perf_2_13_15_Peak/Ref_13_Peak];
Rel_F_2_30 = [Perf_2_6_30_Peak/Ref_6_Peak, Perf_2_8_30_Peak/Ref_8_Peak, ...
    Perf_2_13_30_Peak/Ref_13_Peak];
Rel_F_2_45 = [Perf_2_6_45_Peak/Ref_6_Peak, Perf_2_8_45_Peak/Ref_8_Peak, ...
    Perf_2_13_45_Peak/Ref_13_Peak];

% Perforation 3
Rel_F_3_0 = [Perf_3_6_0_Peak/Ref_6_Peak, Perf_3_8_0_Peak/Ref_8_Peak, ...
    Perf_3_13_0_Peak/Ref_13_Peak];
Rel_F_3_15 = [Perf_3_6_15_Peak/Ref_6_Peak, Perf_3_8_15_Peak/Ref_8_Peak, ...
    Perf_3_13_15_Peak/Ref_13_Peak];
Rel_F_3_30 = [Perf_3_6_30_Peak/Ref_6_Peak, Perf_3_8_30_Peak/Ref_8_Peak, ...
    Perf_3_13_30_Peak/Ref_13_Peak];
Rel_F_3_45 = [Perf_3_6_45_Peak/Ref_6_Peak, Perf_3_8_45_Peak/Ref_8_Peak, ...
   Perf_3_13_45_Peak/Ref_13_Peak];

%% Lines of Best Fit for scaled CFD data
s(1,:) = tanhfit([KC,KC,KC,KC],[Rel_F_1_0,Rel_F_1_15,Rel_F_1_30,Rel_F_1_45]);
s(2,:) = tanhfit([KC,KC,KC,KC],[Rel_F_2_0,Rel_F_2_15,Rel_F_2_30,Rel_F_2_45]);
s(3,:) = tanhfit([KC,KC,KC,KC],[Rel_F_3_0,Rel_F_3_15,Rel_F_3_30,Rel_F_3_45]);

%% Root Mean Square Error (RMSE) in relation to other research
x1 = linspace(0, 20);
bestfit=@(b1,b2,x) b1.*tanh(b2.*x);
RMSE_Ploeg(1) = sqrt(sum((bestfit(1,1.49,x1) - bestfit(s(1,1),s(1,2),x1)).^2)/length(x1));
RMSE_Ploeg(2) = sqrt(sum((bestfit(.93,0.75,x1) - bestfit(s(2,1),s(2,2),x1)) .^2)/length(x1));
RMSE_Ploeg(3) = sqrt(sum((bestfit(0.75,0.52,x1) - bestfit(s(3,1),s(3,2),x1)).^2)/length(x1));

RMSE_Ander(1) = sqrt(sum((bestfit(.94,0.84,x1) - bestfit(s(2,1),s(2,2),x1)).^2)/length(x1));


%% R-Squared Calculation
R_SQ_Ploeg = [R_SQ(bestfit(s(1,1),s(1,2),x1), bestfit(1,1.49,x1));...
    R_SQ(bestfit(s(2,1),s(2,2),x1), bestfit(.93,0.75,x1));...
    R_SQ(bestfit(s(3,1),s(3,2),x1), bestfit(0.75,0.52,x1))];

R_SQ_Anderson = R_SQ(bestfit(.94,0.84,x1), bestfit(s(2,1),s(2,2),x1));



%% Plotting
% Perforation 1
figure(1)
plot(KC,Rel_F_1_0,"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Rel_F_1_15,"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Rel_F_1_30,"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Rel_F_1_45,"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 1*tanh(1.49*x), [1 6], 'g');
fplot(@(x) s(1,1)*tanh(s(1,2)*x), [1 6], 'r');
title("Perforation #1 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim([0.3 1.1])
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Line of Best Fit')
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg (2021)', 'Line of Best Fit')
functionText1 = sprintf('Line of Best Fit: y = %.2f*tanh(%.2f*KC)', s(1,1), s(1,2));
functionText2 = sprintf('RMSE: Ploeg (2021) = %.4f', RMSE_Ploeg(1));
functionText3 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(1));
text(2, .45, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .42, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .38, functionText3, 'FontSize', 10, 'HorizontalAlignment', 'left');

% Perforation 2
figure(2)
plot(KC,Rel_F_2_0,"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Rel_F_2_15,"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Rel_F_2_30,"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Rel_F_2_45,"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 0.94*tanh(0.84*x), [1 6], 'b');
fplot(@(x) 0.93*tanh(0.75*x), [1 6], 'g');
fplot(@(x) s(2,1)*tanh(s(2,2)*x), [1 6], 'r');
title("Perforation #2 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim([0.3 1.1])
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Line of Best Fit')
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Anderson, et al. (2020)', 'Ploeg (2021)', 'Line of Best Fit')
functionText1 = sprintf('Line of Best Fit: y = %.2f*tanh(%.2f*KC)', s(2,1), s(2,2));
functionText2 = sprintf('RMSE: Ploeg (2021) = %.4f', RMSE_Ploeg(2));
functionText3 = sprintf('Anderson, et al. (2020) = %.4f', RMSE_Ander);
functionText4 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(2));
functionText5 = sprintf('Anderson, et al. (2020) = %.4f', R_SQ_Anderson);
text(2, .50, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .47, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2.52, .44, functionText3, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .40, functionText4, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2.25, .37, functionText5, 'FontSize', 10, 'HorizontalAlignment', 'left');

% Perforation 3
figure(3)
plot(KC,Rel_F_3_0,"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Rel_F_3_15,"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Rel_F_3_30,"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Rel_F_3_45,"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 0.75*tanh(0.52*x), [1 6], 'g');
fplot(@(x) s(3,1)*tanh(s(3,2)*x), [1 6], 'r');
title("Perforation #3 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim([0.3 1.1])
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Line of Best Fit')
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg(2021)', 'Line of Best Fit')
functionText1 = sprintf('Line of Best Fit: y = %.2f*tanh(%.2f*KC)', s(3,1), s(3,2));
functionText2 = sprintf('RMSE: Ploeg (2021) = %.4f', RMSE_Ploeg(3));
functionText3 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(3));
text(2, .45, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .42, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .38, functionText3, 'FontSize', 10, 'HorizontalAlignment', 'left');



%% Create table output of all max force values
rownames = {'Reference','Perforation 1', 'Perforation 2', 'Perforation 3'};
columnnames = {'Wave 1', 'Wave 2', 'Wave3'};
Wave1Peaks =[Ref_6_Peak; mean([Perf_1_6_0_Peak,Perf_1_6_15_Peak,Perf_1_6_30_Peak,Perf_1_6_45_Peak]);...
    mean([Perf_2_6_0_Peak,Perf_2_6_15_Peak,Perf_2_6_30_Peak,Perf_2_6_45_Peak]);...
    mean([Perf_3_6_0_Peak,Perf_3_6_15_Peak,Perf_3_6_30_Peak,Perf_3_6_45_Peak])];
Wave2Peaks =[Ref_8_Peak; mean([Perf_1_8_0_Peak,Perf_1_8_15_Peak,Perf_1_8_30_Peak,Perf_1_8_45_Peak]);...
    mean([Perf_2_8_0_Peak,Perf_2_8_15_Peak,Perf_2_8_30_Peak,Perf_2_8_45_Peak]);...
    mean([Perf_3_8_0_Peak,Perf_3_8_15_Peak,Perf_3_8_30_Peak,Perf_3_8_45_Peak])];
Wave3Peaks =[Ref_13_Peak; mean([Perf_1_13_0_Peak,Perf_1_13_15_Peak,Perf_1_13_30_Peak,Perf_1_13_45_Peak]);...
    mean([Perf_2_13_0_Peak,Perf_2_13_15_Peak,Perf_2_13_30_Peak,Perf_2_13_45_Peak]);...
    mean([Perf_3_13_0_Peak,Perf_3_13_15_Peak,Perf_3_13_30_Peak,Perf_3_13_45_Peak])];
T = table(Wave1Peaks, Wave2Peaks, Wave3Peaks, 'RowNames',rownames, 'VariableNames',columnnames)


%% Create table output for rotation based forces
rownames = {'Reference','Perforation 1 - 0 deg', 'Perforation 1 - 15 deg',...
    'Perforation 1 - 30 deg', 'Perforation 1 - 45 deg', 'Perforation 2 - 0 deg',...
    'Perforation 2 - 15 deg', 'Perforation 2 - 30 deg', 'Perforation 2 - 45 deg',...
    'Perforation 3 - 0 deg', 'Perforation 3 - 15 deg', 'Perforation 3 - 30 deg',...
    'Perforation 3 - 45 deg'};
columnnames = {'Wave 1', 'Wave 2', 'Wave3'};
Wave1Peaks =[Ref_6_Peak; Perf_1_6_0_Peak; Perf_1_6_15_Peak; Perf_1_6_30_Peak; Perf_1_6_45_Peak;...
    Perf_2_6_0_Peak; Perf_2_6_15_Peak; Perf_2_6_30_Peak; Perf_2_6_45_Peak;...
    Perf_3_6_0_Peak; Perf_3_6_15_Peak; Perf_3_6_30_Peak; Perf_3_6_45_Peak];
Wave2Peaks =[Ref_8_Peak; Perf_1_8_0_Peak; Perf_1_8_15_Peak; Perf_1_8_30_Peak; Perf_1_8_45_Peak;...
    Perf_2_8_0_Peak; Perf_2_8_15_Peak; Perf_2_8_30_Peak; Perf_2_8_45_Peak;...
    Perf_3_8_0_Peak; Perf_3_8_15_Peak; Perf_3_8_30_Peak; Perf_3_8_45_Peak];
Wave3Peaks =[Ref_13_Peak; Perf_1_13_0_Peak; Perf_1_13_15_Peak; Perf_1_13_30_Peak; Perf_1_13_45_Peak;...
    Perf_2_13_0_Peak; Perf_2_13_15_Peak; Perf_2_13_30_Peak; Perf_2_13_45_Peak;...
    Perf_3_13_0_Peak; Perf_3_13_15_Peak; Perf_3_13_30_Peak; Perf_3_13_45_Peak];
T = table(Wave1Peaks, Wave2Peaks, Wave3Peaks, 'RowNames',rownames, 'VariableNames',columnnames)

%% All perforations
Perf1_Avg_Rel_Force = [mean([Rel_F_1_0(1), Rel_F_1_15(1), Rel_F_1_30(1),Rel_F_1_45(1)]),...
    mean([Rel_F_1_0(2), Rel_F_1_15(2), Rel_F_1_30(2),Rel_F_1_45(2)]),...
    mean([Rel_F_1_0(3), Rel_F_1_15(3), Rel_F_1_30(3),Rel_F_1_45(3)])];
Perf2_Avg_Rel_Force = [mean([Rel_F_2_0(1), Rel_F_2_15(1), Rel_F_2_30(1),Rel_F_2_45(1)]),...
    mean([Rel_F_2_0(2), Rel_F_2_15(2), Rel_F_2_30(2),Rel_F_2_45(2)]),...
    mean([Rel_F_2_0(3), Rel_F_2_15(3), Rel_F_2_30(3),Rel_F_2_45(3)])];
Perf3_Avg_Rel_Force = [mean([Rel_F_3_0(1), Rel_F_3_15(1), Rel_F_3_30(1),Rel_F_3_45(1)]),...
    mean([Rel_F_3_0(2), Rel_F_3_15(2), Rel_F_3_30(2),Rel_F_3_45(2)]),...
    mean([Rel_F_3_0(3), Rel_F_3_15(3), Rel_F_3_30(3),Rel_F_3_45(3)])];
figure(4)
plot(KC,Perf1_Avg_Rel_Force(:),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf2_Avg_Rel_Force(:),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf3_Avg_Rel_Force(:),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
fplot(@(x) s(1,1)*tanh(s(1,2)*x), [1 6], 'r');
fplot(@(x) s(2,1)*tanh(s(2,2)*x), [1 6], 'g');
fplot(@(x) s(3,1)*tanh(s(3,2)*x), [1 6], 'b');
ylim([0.3 1.1])
title("Relative Force vs. KC number - CFD Results")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
legend('Perforation 1', 'Perforation 2', 'Perforation 3', 'Perforation 1 Line of Best Fit', 'Perforation 2 Line of Best Fit','Perforation 3 Line of Best Fit')

%% Plotting all perforations with rotations
% Perforation 1
figure()
plot([KC, KC, KC],[Rel_F_1_0, Rel_F_2_0, Rel_F_3_0],"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot([KC, KC, KC],[Rel_F_1_15, Rel_F_2_15, Rel_F_3_15],"Color",[0 1 0],"Marker","o","LineStyle","none")
plot([KC, KC, KC],[Rel_F_1_30, Rel_F_2_30, Rel_F_3_30],"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot([KC, KC, KC],[Rel_F_1_45, Rel_F_2_45, Rel_F_3_45],"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) s(1,1)*tanh(s(1,2)*x), [1 6], 'r');
fplot(@(x) s(2,1)*tanh(s(2,2)*x), [1 6], 'g');
fplot(@(x) s(3,1)*tanh(s(3,2)*x), [1 6], 'b');

title("Relative Force vs. KC number - CFD Results with Rotation")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim([0.3 1.1])
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', ...
    'Perforation 1 Line of Best Fit', 'Perforation 2 Line of Best Fit',...
    'Perforation 3 Line of Best Fit')


