%% Wave Tank Parameters Scaled Monopile Data Analysis
% For the data analysis of forces on scaled monopiles with and without
% perforations. This script will output 3 graphs for each of the
% perforation geometries relating the force on the perforation/force on the
% reference vs. KC number. Each graph contains the trends for 1 wave
% approach angle (0) at 3 wave test states (1, 2, 3 from wave tank data) and
% the applicable trend from either Anderson, et al. (2020) and Ploeg (2021)

%%RUN AFTER RUNNING WAVE TANK ANALYSIS 2

% clc
% clear
% close all

addpath("Import\");

%% Variable Setup
KC = [1.76474,4.739643,3.001749];
datafilepath = "H:\Thesis\Scaled Simulations\Wave Tank Parameters\";
datafilepath2 = "H:\Thesis\Scaled Simulations\Anderson\";

%% Data Import
% Reference (3 total)
Ref_1 = table2array(import01(datafilepath + "Reference_1.txt", 9, inf));
Ref_2 =  table2array(importPressAndForce(datafilepath + "Reference_2.txt", 9, inf));
Ref_3 =  table2array(import01(datafilepath + "Reference_3.txt", 9, inf));

% Perforation 1 (12 total)
Perf_1_1_0 =  table2array(import01(datafilepath + "Perforation_1_1_0.txt", 9, inf));
Perf_1_2_0 = table2array(import01(datafilepath + "Perforation_1_2_0.txt", 9, inf));
Perf_1_3_0 = table2array(import01(datafilepath + "Perforation_1_3_0.txt", 9, inf));

% Perforation 2 (12 total)
Perf_2_1_0 = table2array(import01(datafilepath + "Perforation_2_1_0.txt", 9, inf));
Perf_2_2_0 = table2array(import01(datafilepath + "Perforation_2_2_0.txt", 9, inf));
Perf_2_3_0 = table2array(import01(datafilepath + "Perforation_2_3_0.txt", 9, inf));

% Perforation 3 (12 total).
Perf_3_1_0 = table2array(import01(datafilepath + "Perforation_3_1_0.txt", 9, inf));
Perf_3_2_0 = table2array(import01(datafilepath + "Perforation_3_2_0.txt", 9, inf));
Perf_3_3_0 = table2array(import01(datafilepath + "Perforation_3_3_0.txt", 9, inf));

Ref_2_Ander =  table2array(importPressAndForce(datafilepath2 + "Reference_2-3.txt", 9, inf));
Perf_2_2_Ander =  table2array(importPressAndForce(datafilepath2 + "Perforation_2_2_0.txt", 9, inf));


%% Average of Peak Total Force
% Reference
Ref_1_Peak = mean(maxk(Ref_1(4:150,3),6,1));
Ref_2_Peak = mean(maxk(Ref_2(3:510,3),6,1));
Ref_3_Peak = mean(maxk(Ref_3(3:end,3),2,1));

% Perforation 1
Perf_1_1_0_Peak = mean(maxk(Perf_1_1_0(3:200,3),6,1));
Perf_1_2_0_Peak = mean(maxk(Perf_1_2_0(3:510,3),6,1));
Perf_1_3_0_Peak = mean(maxk(Perf_1_3_0(3:530,3),6,1));

% Perfotation 2
Perf_2_1_0_Peak = mean(maxk(Perf_2_1_0(3:200,3),6,1));
Perf_2_2_0_Peak = mean(maxk(Perf_2_2_0(3:510,3),6,1));
Perf_2_3_0_Peak = mean(maxk(Perf_2_3_0(3:530,3),6,1));

% Perforation 3
Perf_3_1_0_Peak = mean(maxk(Perf_3_1_0(3:200,3),6,1));
Perf_3_2_0_Peak = mean(maxk(Perf_3_2_0(3:510,3),6,1));
Perf_3_3_0_Peak = mean(maxk(Perf_3_3_0(3:530,3),6,1));

% Anderson
Ref_2_Ander_Peak = mean(maxk(Ref_2_Ander(3:510,3),6,1));
Perf_2_2_0_Ander_Peak = mean(maxk(Perf_2_2_Ander(3:510,3),6,1));

%% Relation between Perforation Force and Total Force
% Perforation 1
Rel_F_1_0 = [Perf_1_1_0_Peak/Ref_1_Peak, Perf_1_2_0_Peak/Ref_2_Peak, ...
    Perf_1_3_0_Peak/Ref_3_Peak];

% Perforation 2
Rel_F_2_0 = [Perf_2_1_0_Peak/Ref_1_Peak, Perf_2_2_0_Peak/Ref_2_Peak, ...
    Perf_2_3_0_Peak/Ref_3_Peak];

% Perforation 3
Rel_F_3_0 = [Perf_3_1_0_Peak/Ref_1_Peak, Perf_3_2_0_Peak/Ref_2_Peak, ...
    Perf_3_3_0_Peak/Ref_3_Peak];

%% Lines of Best Fit for scaled CFD data
s(1,:) = tanhfit([KC(1) KC(3) KC(2)],[Rel_F_1_0(1) Rel_F_1_0(3) Rel_F_1_0(2)]);
s(2,:) = tanhfit([KC(1) KC(3) KC(2)],[Rel_F_2_0(1) Rel_F_2_0(3) Rel_F_2_0(2)]);
s(3,:) = tanhfit([KC(1) KC(3) KC(2)],[Rel_F_3_0(1) Rel_F_3_0(3) Rel_F_3_0(2)]);

%% Root Mean Square Error (RMSE) in relation to previous CFD analysis
x1 = linspace(0, 20);
bestfit=@(b1,b2,x) b1.*tanh(b2.*x);
RMSE_WT(1) = sqrt(sum((bestfit(s2(1,1),s2(1,2),x1) - bestfit(s(1,1),s(1,2),x1)).^2)/length(x1));
RMSE_WT(2) = sqrt(sum((bestfit(s2(2,1),s2(2,2),x1) - bestfit(s(2,1),s(2,2),x1)) .^2)/length(x1));
RMSE_WT(3) = sqrt(sum((bestfit(s2(3,1),s2(3,2),x1) - bestfit(s(3,1),s(3,2),x1)).^2)/length(x1));

%% R-Squared Calculation
R_SQ_WT(1) = R_SQ(bestfit(s(1,1),s(1,2),x1), bestfit(s2(1,1),s2(1,2),x1));
R_SQ_WT(2) = R_SQ(bestfit(s(2,1),s(2,2),x1), bestfit(s2(2,1),s2(2,2),x1));
R_SQ_WT(3) = R_SQ(bestfit(s(3,1),s(3,2),x1), bestfit(s2(3,1),s2(3,2),x1));


% %% RMSE to Previous CFD Analysis
% RMSE_WT(1) = sqrt(sum([mean((Perf1_Rel_Force(1,:))-Rel_F_1_0(1))^2, ...
%     (mean(Perf1_Rel_Force(2,:))-Rel_F_1_0(2))^2, ...
%     (mean(Perf1_Rel_Force(3,:))-Rel_F_1_0(3))^2])/3);
% RMSE_WT(2) = sqrt(sum([mean((Perf2_Rel_Force(1,:))-Rel_F_2_0(1))^2, ...
%     (mean(Perf2_Rel_Force(2,:))-Rel_F_2_0(2))^2, ...
%     (mean(Perf2_Rel_Force(3,:))-Rel_F_2_0(3))^2])/3);
% RMSE_WT(3) = sqrt(sum([mean((Perf3_Rel_Force(1,:))-Rel_F_3_0(1))^2, ...
%     (mean(Perf3_Rel_Force(2,:))-Rel_F_3_0(2))^2, ...
%     (mean(Perf3_Rel_Force(3,:))-Rel_F_3_0(3))^2])/3);
% 
% %% R-Squared Calculation
% R_SQ_WT(1) = R_SQ([mean(Perf1_Rel_Force(1,:)), mean(Perf1_Rel_Force(2,:)), mean(Perf1_Rel_Force(3,:))],...
%     [Rel_F_1_0(1), Rel_F_1_0(2), Rel_F_1_0(3)]);
% R_SQ_WT(2) = R_SQ([mean(Perf2_Rel_Force(1,:)), mean(Perf2_Rel_Force(2,:)), mean(Perf2_Rel_Force(3,:))],...
%     [Rel_F_3_0(1), Rel_F_2_0(2), Rel_F_2_0(3)]);
% R_SQ_WT(3) = R_SQ([mean(Perf3_Rel_Force(1,:)), mean(Perf3_Rel_Force(2,:)), mean(Perf3_Rel_Force(3,:))],...
%     [Rel_F_3_0(1), Rel_F_3_0(2), Rel_F_3_0(3)]);



%% Plotting
% Perforation 1
figure(1)
plot(KC,Perf1_Rel_Force(:,1),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf1_Rel_Force(:,2),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf1_Rel_Force(:,3),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Perf1_Rel_Force(:,4),"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 1*tanh(1.49*x), [1 6], 'g');
fplot(@(x) 0.99*tanh(0.92*x), [1 6], 'r');
plot(KC,Rel_F_1_0,"Color",[1 0 0],"Marker","*","LineStyle","none")
title("Perforation #1 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim ([0.3 1.1])
legend('0 Degrees','15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg (2021)', 'CFD Line of Best Fit', 'CFD Retest')
functionText1 = sprintf('RMSE = %.4f', RMSE_WT(1));
functionText2 = sprintf('R^2 = %.4f', R_SQ_WT(1));
text(1.5, .4, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(1.5, .36, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');

% Perforation 2
figure(2)
plot(KC,Perf2_Rel_Force(:,1),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf2_Rel_Force(:,2),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf2_Rel_Force(:,3),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Perf2_Rel_Force(:,4),"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 0.94*tanh(0.84*x), [1 6], 'b');
fplot(@(x) 0.93*tanh(0.75*x), [1 6], 'g');
fplot(@(x) 0.91*tanh(0.67*x), [1 6], 'r');
plot(KC,Rel_F_2_0,"Color",[1 0 0],"Marker","*","LineStyle","none")
title("Perforation #2 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim ([0.3 1.1])
legend('0 Degrees','15 Degrees', '30 Degrees', '45 Degrees', 'Anderson, et al.(2020)', 'Ploeg (2021)', 'CFD Line of Best Fit', 'CFD Retest')
functionText1 = sprintf('RMSE = %.4f', RMSE_WT(2));
functionText2 = sprintf('R^2 = %.4f', R_SQ_WT(2));
text(1.5, .4, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(1.5, .36, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');

% Perforation 3
figure(3)
plot(KC,Perf3_Rel_Force(:,1),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf3_Rel_Force(:,2),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf3_Rel_Force(:,3),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Perf3_Rel_Force(:,4),"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 0.75*tanh(0.52*x), [1 6], 'g');
fplot(@(x) 0.68*tanh(0.62*x), [1 6], 'r');
plot(KC,Rel_F_3_0,"Color",[1 0 0],"Marker","*","LineStyle","none")
title("Perforation #3 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
ylim ([0.3 1.1])
legend('0 Degrees','15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg (2021)', 'CFD Line of Best Fit', 'CFD Retest')
functionText1 = sprintf('RMSE = %.4f', RMSE_WT(3));
functionText2 = sprintf('R^2 = %.4f', R_SQ_WT(3));
text(2, .4, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2, .36, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');



%% Create table output of all max force values
rownames = {'Reference','Perforation 1', 'Perforation 2', 'Perforation 3'};
columnnames = {'Wave 1', 'Wave 2', 'Wave3'};
Wave1Peaks =[Ref_1_Peak; Perf_1_1_0_Peak; Perf_2_1_0_Peak; Perf_3_1_0_Peak];
Wave2Peaks =[Ref_2_Peak; Perf_1_2_0_Peak; Perf_2_2_0_Peak; Perf_3_2_0_Peak];
Wave3Peaks =[Ref_3_Peak; Perf_1_3_0_Peak; Perf_2_3_0_Peak; Perf_3_3_0_Peak];
T = table(Wave1Peaks, Wave2Peaks, Wave3Peaks, 'RowNames',rownames, 'VariableNames',columnnames)