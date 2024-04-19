%% CFD Wave Param 
% CFD Wave Parameter Check for WT Param Input

clc
clear
close all 

addpath("Import\");

%% Import Wave Elevations
Ref_1_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Reference_1_Eta.txt", [10, Inf]));
Ref_2_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Reference_2_Eta.txt", [10, Inf]));
Ref_3_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Reference_3_Eta.txt", [10, Inf]));

Perf1_1_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_1_1_0_Eta.txt", [10, Inf]));
Perf1_2_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_1_2_0_Eta_Mesh1.txt", [10, Inf]));
Perf1_3_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_1_3_0_Eta.txt", [10, Inf]));

Perf2_1_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_2_1_0_Eta.txt", [10, Inf]));
Perf2_2_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_2_2_0_Eta_Mesh1.txt", [10, Inf]));
Perf2_3_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_2_3_0_Eta.txt", [10, Inf]));

Perf3_1_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_3_1_0_Eta.txt", [10, Inf]));
Perf3_2_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_3_2_0_Eta_Mesh1.txt", [10, Inf]));
Perf3_3_0_Eta = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_3_3_0_Eta.txt", [10, Inf]));

Slice = table2array(importEta("H:\Thesis\Scaled Simulations\Wave Tank Parameters\Slice-Ander.txt", [10, Inf]));

Ref_2_Ander = table2array(importEta("H:\Thesis\Scaled Simulations\Anderson\Reference_2-3_Eta.txt", [10, Inf]));
Perf_2_2_Ander = table2array(importEta("H:\Thesis\Scaled Simulations\Anderson\Perforation_2_2_0_Eta.txt", [10, Inf]));

%% Pass through zeroup analysis
for n = 1:12
    Avg(n).Wave= n-3*(floor(n/3.1));
end
[~,~,~, Avg(1).H, Avg(1).T] = ZeroUp((Ref_1_Eta(:,2)-.4375),Ref_1_Eta(:,1));
[~, Avg(1).L] = WaveSolver(Avg(1).T,.4375,Avg(1).H);
[~,~,~, Avg(2).H, Avg(2).T] = ZeroUp((Ref_2_Eta(1:55,2)-.4375),Ref_2_Eta(1:55,1));
[~, Avg(2).L] = WaveSolver(Avg(2).T,.4375,Avg(2).H);
[~,~,~, Avg(3).H, Avg(3).T] = ZeroUp((Ref_3_Eta(1:53,2)-.4375),Ref_3_Eta(1:53,1));
[~, Avg(3).L] = WaveSolver(Avg(3).T,.4375,Avg(3).H);

[~,~,~, Avg(4).H, Avg(4).T] = ZeroUp((Perf1_1_0_Eta(:,2)-.4375),Perf1_1_0_Eta(:,1));
[~, Avg(4).L] = WaveSolver(Avg(4).T,.4375,Avg(4).H);
[~,~,~, Avg(5).H, Avg(5).T] = ZeroUp((Perf1_2_0_Eta(1:55,2)-.4375),Perf1_2_0_Eta(1:55,1));
[~, Avg(5).L] = WaveSolver(Avg(5).T,.4375,Avg(5).H);
[~,~,~, Avg(6).H, Avg(6).T] = ZeroUp((Perf1_3_0_Eta(1:53,2)-.4375),Perf1_3_0_Eta(1:53,1));
[~, Avg(6).L] = WaveSolver(Avg(6).T,.4375,Avg(6).H);

[~,~,~, Avg(7).H, Avg(7).T] = ZeroUp((Perf2_1_0_Eta(:,2)-.4375),Perf2_1_0_Eta(:,1));
[~, Avg(7).L] = WaveSolver(Avg(7).T,.4375,Avg(7).H);
[~,~,~, Avg(8).H, Avg(8).T] = ZeroUp((Perf2_2_0_Eta(1:55,2)-.4375),Perf2_2_0_Eta(1:55,1));
[~, Avg(8).L] = WaveSolver(Avg(8).T,.4375,Avg(8).H);
[~,~,~, Avg(9).H, Avg(9).T] = ZeroUp((Perf2_3_0_Eta(1:53,2)-.4375),Perf2_3_0_Eta(1:53,1));
[~, Avg(9).L] = WaveSolver(Avg(9).T,.4375,Avg(9).H);

[~,~,~, Avg(10).H, Avg(10).T] = ZeroUp((Perf3_1_0_Eta(:,2)-.4375),Perf3_1_0_Eta(:,1));
[~, Avg(10).L] = WaveSolver(Avg(10).T,.4375,Avg(10).H);
[~,~,~, Avg(11).H, Avg(11).T] = ZeroUp((Perf3_2_0_Eta(1:75,2)-.4375),Perf3_2_0_Eta(1:75,1));
[~, Avg(11).L] = WaveSolver(Avg(11).T,.4375,Avg(11).H);
[~,~,~, Avg(12).H, Avg(12).T] = ZeroUp((Perf3_3_0_Eta(1:53,2)-.4375),Perf3_3_0_Eta(1:53,1));
[~, Avg(12).L] = WaveSolver(Avg(12).T,.4375,Avg(12).H);

[~,~,~, AvgH, AvgT] = ZeroUp((Slice(1:60,2)-.4375),Slice(1:60,1));

[~,~,~, AvgH2, AvgT2] = ZeroUp((Ref_2_Ander(1:60,2)-.4375),Ref_2_Ander(1:60,1));
[~,~,~, AvgH3, AvgT3] = ZeroUp((Perf_2_2_Ander(1:60,2)-.4375),Perf_2_2_Ander(1:60,1));

%% Plot Like Waves with Like
figure()
plot(Ref_2_Ander(1:60,1), Ref_2_Ander(1:60,2)-.4375)
hold on
