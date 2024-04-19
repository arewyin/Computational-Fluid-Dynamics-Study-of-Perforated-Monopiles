%% Wave Tank Analysis 3
% Pulls data from wave gauges and load cell for analysis with formatting to
% reduce the size of the program line wise. This version uses structures to
% hold the data for post processing instead of individual variables for
% everything
% 
clc
clear 
close all
addpath("Import\");

%%
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontSize', 12) 
set(0,'defaultAxesFontName','Arial')

%% Set Working DRIVE
DRIVE = "h:\";
SCALED_TESTS_ROOT_DIR = DRIVE + "Thesis\Wave Tank Tests\Scaled ReTest\";

%% Load in Paths
%Reference Monopile
Wave_Gauge_Info_Ref = importRefPaths(SCALED_TESTS_ROOT_DIR +...
    "Wave_Gauge_Paths_Ref.csv", [2, Inf]);
Wave_Gauge_Paths_Ref = SCALED_TESTS_ROOT_DIR + Wave_Gauge_Info_Ref.Path;
Wave_Gauge_Limits_Ref = [Wave_Gauge_Info_Ref.Start, Wave_Gauge_Info_Ref.End,...
    Wave_Gauge_Info_Ref.Wave, Wave_Gauge_Info_Ref.Run];
Load_Cell_Info_Ref = importRefLoadPath(SCALED_TESTS_ROOT_DIR+...
    "Wave_Load_Cell_Paths_Ref.csv", [2, Inf]);
Load_Cell_Paths_Ref = SCALED_TESTS_ROOT_DIR + Load_Cell_Info_Ref.Path;
Load_Cell_Limits_Ref = [Load_Cell_Info_Ref.Start, Load_Cell_Info_Ref.End,...
    Load_Cell_Info_Ref.DataStart, Load_Cell_Info_Ref.DataEnd, ...
    Load_Cell_Info_Ref.Wave,Load_Cell_Info_Ref.Run];

%Perforation #2
Wave_Gauge_Paths_Perf2 = SCALED_TESTS_ROOT_DIR + "20240315-Perf2-2-0-7";
Wave_Gauge_Limits_Perf2 = [775, 1127, 2, 1];
Load_Cell_Paths_Perf2 = SCALED_TESTS_ROOT_DIR + "20240315-Perf2-2-0-7-Force";
Load_Cell_Limits_Perf2 = [5, 50, 86, 117, 2,1];

%% Ref Wave Gauge Data Load In
%Initialize locations table
sz = [1,4];
Vartypes = ["double", "double", "double", "double"];
Varnames = ["Wave #", "Run #", "Gauge #", "Location"];
Ref_Loc = table('Size', sz, 'VariableTypes',Vartypes,'VariableNames',Varnames);

%Initialize location
n = 1;

%Loop through paths and assign data to a nonscalar stucture
for i = 1:length(Wave_Gauge_Paths_Ref)
    %Initialize Variables
    Path = Wave_Gauge_Paths_Ref(i,1);
    Path1 = Load_Cell_Paths_Ref(i);
    SaveFile = sprintf('Ref-Wave %d-Run %d',Wave_Gauge_Limits_Ref(i,3),...
        Wave_Gauge_Limits_Ref(i,4));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] = ...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Ref(i,1),...
        Wave_Gauge_Limits_Ref(i,2), SaveFile, SavePath);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(1).T), 0.4375, Wave(1).H);      %Determine Wave Length For Individual Waves
        % Build Nonscalar Structure
        Ref(n).WaveNum = Wave_Gauge_Limits_Ref(i,3);
        Ref(n).RunNum = Wave_Gauge_Limits_Ref(i,4);
        Ref(n).GaugeNum = k;
        Ref(n).WaveTime = WavesTime(:);
        Ref(n).WaveLimits = Wave(k).Waves;
        Ref(n).WaveEta = Waves(:,k);
        Ref(n).WaveHeight = AvgH(k);
        Ref(n).WaveHm0 = Hm0(k);
        Ref(n).WavePeriod = AvgT(k);
        Ref(n).WaveLength = L;
        Ref(n).Ursells = AvgH(k)/100*L^2/(0.4375^3);
        Ref_Loc(n,:) = {Ref(n).WaveNum, Ref(n).RunNum, Ref(n).GaugeNum, n};
        n = n+1;
    end
    

    % Determine max force
    [ForceSeries, AveragedForces, MaxForce] = ...
    Load_Cell_Force(Path1, Load_Cell_Limits_Ref(i,1), ...
    Load_Cell_Limits_Ref(i,2), Load_Cell_Limits_Ref(i,3), ...
    Load_Cell_Limits_Ref(i,4),  L2, Wave(3).H, SaveFile, SavePath);


    % Create structure to hold force time series forces averaged and
    % max force for a given wave and run.
    Ref_Force(i).WaveNum = Load_Cell_Limits_Ref(i,5);
    Ref_Force(i).RunNum = Load_Cell_Limits_Ref(i,6);
    Ref_Force(i).ForceSeries = ForceSeries;
    Ref_Force(i).AveragedForces = AveragedForces;
    Ref_Force(i).MaxForce = MaxForce;

    clear WavesTime
    clear Waves
    clear ForceSeries
    clear AveragedForces
    clear MaxForce
end


%% Perf2 Wave Gauge Data Load In
%Initialize location
n = 1;

%Loop through paths and assign data to a nonscalar stucture
for i = 1:length(Wave_Gauge_Paths_Perf2)
    %Initialize Variables
    Path = Wave_Gauge_Paths_Perf2(i,1);
    Path1 = Load_Cell_Paths_Perf2(i);
    SaveFile = sprintf('Perf2-Wave %d-Run %d',Wave_Gauge_Limits_Perf2(i,3),...
        Wave_Gauge_Limits_Perf2(i,4));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] = ...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Perf2(i,1),...
        Wave_Gauge_Limits_Perf2(i,2), SaveFile, SavePath, 2);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(1).T), 0.4375, Wave(1).H);      %Determine Wave Length For Individual Waves
        % Build Nonscalar Structure
        Perf2(n).WaveNum = Wave_Gauge_Limits_Perf2(i,3);
        Perf2(n).RunNum = Wave_Gauge_Limits_Perf2(i,4);
        Perf2(n).GaugeNum = k;
        Perf2(n).WaveTime = WavesTime(:);
        Perf2(n).WaveLimits = Wave(k).Waves;
        Perf2(n).WaveEta = Waves(:,k);
        Perf2(n).WaveHeight = AvgH(k);
        Perf2(n).WaveHm0 = Hm0(k);
        Perf2(n).WavePeriod = AvgT(k);
        Perf2(n).WaveLength = L;
        Perf2(n).Ursells = AvgH(k)/100*L^2/(0.4375^3);
        n = n+1;
    end
    

    % Determine max force
    [ForceSeries, AveragedForces, MaxForce] = ...
    Load_Cell_Force(Path1, Load_Cell_Limits_Perf2(i,1), ...
    Load_Cell_Limits_Perf2(i,2), Load_Cell_Limits_Perf2(i,3), ...
    Load_Cell_Limits_Perf2(i,4),  L2, Wave(3).H, SaveFile, SavePath);


    % Create structure to hold force time series forces averaged and
    % max force for a given wave and run.
    Perf2_Force(i).WaveNum = Load_Cell_Limits_Ref(i,5);
    Perf2_Force(i).RunNum = Load_Cell_Limits_Ref(i,6);
    Perf2_Force(i).ForceSeries = ForceSeries;
    Perf2_Force(i).AveragedForces = AveragedForces;
    Perf2_Force(i).MaxForce = MaxForce;

    clear WavesTime
    clear Waves
    clear ForceSeries
    clear AveragedForces
    clear MaxForce
end




%% Create table output of all max force values
rownames = {'Max Force','Wave Height', 'Wave Period', 'Wave Length'};
columnnames = {'Run 1', 'Run 2', 'Run 3'};
Run1 =[Ref_Force(1).MaxForce; Ref(1).WaveHeight; Ref(1).WavePeriod; Ref(1).WaveLength];
Run2 =[Ref_Force(2).MaxForce; Ref(4).WaveHeight; Ref(4).WavePeriod; Ref(4).WaveLength];
Run3 =[Ref_Force(3).MaxForce; Ref(7).WaveHeight; Ref(7).WavePeriod; Ref(7).WaveLength];
T1 = table(Run1, Run2, Run3, 'RowNames',rownames, 'VariableNames',columnnames)


%% Create table output of all max force values
rownames = {'Max Force','Wave Height', 'Wave Period', 'Wave Length'};
columnnames = {'Perforation #1'};
Run1 =[Perf2_Force(1).MaxForce; Perf2(1).WaveHeight; Perf2(1).WavePeriod; Perf2(1).WaveLength];
T1 = table(Run1, 'RowNames',rownames, 'VariableNames',columnnames)