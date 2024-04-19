%% Wave Tank Analysis 2
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
KC = [1.75,4.73,3.0];

%% Set Working DRIVE
DRIVE = "F:\";
SCALED_TESTS_ROOT_DIR = DRIVE + "Thesis\Wave Tank Tests\Scaled Tests\";

%% Load in Paths
%Reference Monopile
Wave_Gauge_Info_Ref = importRefPaths(SCALED_TESTS_ROOT_DIR +...
    "Wave_Gauge_Paths_Ref.csv", [2, Inf]);
Wave_Gauge_Paths_Ref = Wave_Gauge_Info_Ref.Path;
Wave_Gauge_Limits_Ref = [Wave_Gauge_Info_Ref.Start, Wave_Gauge_Info_Ref.End,...
    Wave_Gauge_Info_Ref.Wave, Wave_Gauge_Info_Ref.Run];
Load_Cell_Info_Ref = importRefLoadPath(SCALED_TESTS_ROOT_DIR+...
    "Wave_Load_Cell_Paths_Ref.csv", [2, Inf]);
Load_Cell_Paths_Ref = Load_Cell_Info_Ref.Path;
Load_Cell_Limits_Ref = [Load_Cell_Info_Ref.Start, Load_Cell_Info_Ref.End,...
    Load_Cell_Info_Ref.DataStart, Load_Cell_Info_Ref.DataEnd, ...
    Load_Cell_Info_Ref.Wave,Load_Cell_Info_Ref.Run];

%Perfoation #1
Wave_Gauge_Info_Perf1 = importPerfPath(SCALED_TESTS_ROOT_DIR +...
    "Wave_Gauge_Paths_Perf1.csv", [2, Inf]);
Wave_Gauge_Paths_Perf1 = Wave_Gauge_Info_Perf1.Path;
Wave_Gauge_Limits_Perf1 = [Wave_Gauge_Info_Perf1.Start, Wave_Gauge_Info_Perf1.End,...
    Wave_Gauge_Info_Perf1.Wave, Wave_Gauge_Info_Perf1.Rotation, Wave_Gauge_Info_Perf1.Run];
Load_Cell_Info_Perf1 = importPerfLoadPath(SCALED_TESTS_ROOT_DIR +...
    "Wave_Load_Cell_Paths_Perf1.csv", [2, Inf]);
Load_Cell_Paths_Perf1 = Load_Cell_Info_Perf1.Path;
Load_Cell_Limits_Perf1 = [Load_Cell_Info_Perf1.Start, Load_Cell_Info_Perf1.End,...
    Load_Cell_Info_Perf1.DataStart, Load_Cell_Info_Perf1.DataEnd, ...
    Load_Cell_Info_Perf1.Wave, Load_Cell_Info_Perf1.Rotation,...
    Load_Cell_Info_Perf1.Run];

%Perfoation #2
Wave_Gauge_Info_Perf2 = importPerfPath(SCALED_TESTS_ROOT_DIR+...
    "Wave_Gauge_Paths_Perf2.csv", [2, Inf]);
Wave_Gauge_Paths_Perf2 = Wave_Gauge_Info_Perf2.Path;
Wave_Gauge_Limits_Perf2 = [Wave_Gauge_Info_Perf2.Start, Wave_Gauge_Info_Perf2.End,...
    Wave_Gauge_Info_Perf2.Wave, Wave_Gauge_Info_Perf2.Rotation, Wave_Gauge_Info_Perf2.Run];
Load_Cell_Info_Perf2 = importPerfLoadPath(SCALED_TESTS_ROOT_DIR +...
    "Wave_Load_Cell_Paths_Perf2.csv", [2, Inf]);
Load_Cell_Paths_Perf2 = Load_Cell_Info_Perf2.Path;
Load_Cell_Limits_Perf2 = [Load_Cell_Info_Perf2.Start, Load_Cell_Info_Perf2.End,...
    Load_Cell_Info_Perf2.DataStart, Load_Cell_Info_Perf2.DataEnd, ...
    Load_Cell_Info_Perf2.Wave, Load_Cell_Info_Perf2.Rotation,...
    Load_Cell_Info_Perf2.Run];

%Perfoation #3
Wave_Gauge_Info_Perf3 = importPerfPath(SCALED_TESTS_ROOT_DIR+...
    "Wave_Gauge_Paths_Perf3.csv", [2, Inf]);
Wave_Gauge_Paths_Perf3 = Wave_Gauge_Info_Perf3.Path;
Wave_Gauge_Limits_Perf3 = [Wave_Gauge_Info_Perf3.Start, Wave_Gauge_Info_Perf3.End,...
    Wave_Gauge_Info_Perf3.Wave, Wave_Gauge_Info_Perf3.Rotation, Wave_Gauge_Info_Perf3.Run];
Load_Cell_Info_Perf3 = importPerfLoadPath(SCALED_TESTS_ROOT_DIR +...
    "Wave_Load_Cell_Paths_Perf3.csv", [2, Inf]);
Load_Cell_Paths_Perf3 = Load_Cell_Info_Perf3.Path;
Load_Cell_Limits_Perf3 = [Load_Cell_Info_Perf3.Start, Load_Cell_Info_Perf3.End,...
    Load_Cell_Info_Perf3.DataStart, Load_Cell_Info_Perf3.DataEnd, ...
    Load_Cell_Info_Perf3.Wave, Load_Cell_Info_Perf3.Rotation,...
    Load_Cell_Info_Perf3.Run];

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
    Path = strcat(DRIVE, Wave_Gauge_Paths_Ref(i,1));
    Path1 = strcat(DRIVE, Load_Cell_Paths_Ref(i));
    SaveFile = sprintf('Ref-Wave %d-Run %d',Wave_Gauge_Limits_Ref(i,3),...
        Wave_Gauge_Limits_Ref(i,4));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] = ...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Ref(i,1),...
        Wave_Gauge_Limits_Ref(i,2), SaveFile, SavePath);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(3).T), 0.4375, Wave(3).H);      %Determine Wave Length For Individual Waves
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


%% Perf1 Wave Gauge Data Load In
%Initialize locations table
sz = [1,5];
Vartypes = ["double", "double", "double", "double", "double"];
Varnames = ["Wave #", "Run #", "Rotation", "Gauge #", "Location"];
Perf1_Loc = table('Size', sz, 'VariableTypes',Vartypes,'VariableNames',Varnames);

%Initialize location
n = 1;

%Loop through paths and assign data to a nonscalar stucture
for i = 1:length(Wave_Gauge_Paths_Perf1)
    %Initialize Variables
    Path = strcat(DRIVE, Wave_Gauge_Paths_Perf1(i,1));
    Path1 = strcat(DRIVE, Load_Cell_Paths_Perf1(i));
    SaveFile = sprintf('Perf 1-Wave %d-%d Degrees-Run %d',Wave_Gauge_Limits_Perf1(i,3),...
        Wave_Gauge_Limits_Perf1(i,4), Wave_Gauge_Limits_Perf1(i,5));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] = ...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Perf1(i,1),...
        Wave_Gauge_Limits_Perf1(i,2), SaveFile, SavePath);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(3).T), 0.4375, Wave(3).H);      %Determine Wave Length For Individual Waves
        % Build Nonscalar Structure
        Perf1(n).WaveNum = Wave_Gauge_Limits_Perf1(i,3);
        Perf1(n).Rotation = Wave_Gauge_Limits_Perf1(i,4);
        Perf1(n).RunNum = Wave_Gauge_Limits_Perf1(i,5);
        Perf1(n).GaugeNum = k;
        Perf1(n).WaveTime = WavesTime(:);
        Perf1(n).WaveLimits = Wave(k).Waves;
        Perf1(n).WaveEta = Waves(:,k);
        Perf1(n).WaveHeight = AvgH(k);
        Perf1(n).WaveHm0 = Hm0(k);
        Perf1(n).WavePeriod = AvgT(k);
        Perf1(n).WaveLength = L;
        Perf1(n).Ursells = AvgH(k)/100*L^2/(0.4375^3);
        Perf1_Loc(n,:) = {Perf1(n).WaveNum, Perf1(n).Rotation,...
            Perf1(n).RunNum, Perf1(n).GaugeNum, n};
        n = n+1;
    end

    % Determine max force
    [ForceSeries, AveragedForces, MaxForce] = ...
    Load_Cell_Force(Path1, Load_Cell_Limits_Perf1(i,1), ...
    Load_Cell_Limits_Perf1(i,2), Load_Cell_Limits_Perf1(i,3), ...
    Load_Cell_Limits_Perf1(i,4), L2, Wave(3).H, SaveFile, SavePath);


    % Create structure to hold force time series forces averaged and
    % max force for a given wave and run.
    Perf1_Force(i).WaveNum = Load_Cell_Limits_Perf1(i,5);
    Perf1_Force(i).RunNum = Load_Cell_Limits_Perf1(i,7);
    Perf1_Force(i).Rotation = Load_Cell_Limits_Perf1(i,6);
    Perf1_Force(i).ForceSeries = ForceSeries;
    Perf1_Force(i).AveragedForces = AveragedForces;
    Perf1_Force(i).MaxForce = MaxForce;

    clear WavesTime
    clear Waves
    clear ForceSeries
    clear AveragedForces
    clear MaxForce
end

%% Perf2 Wave Gauge Data Load In
%Initialize locations table
sz = [1,5];
Vartypes = ["double", "double", "double", "double", "double"];
Varnames = ["Wave #", "Run #", "Rotation", "Gauge #", "Location"];
Perf2_Loc = table('Size', sz, 'VariableTypes',Vartypes,'VariableNames',Varnames);

%Initialize location
n = 1;

%Loop through paths and assign data to a nonscalar stucture
for i = 1:length(Wave_Gauge_Paths_Perf2)
    %Initialize Variables
    Path = strcat(DRIVE,Wave_Gauge_Paths_Perf2(i,1));
    Path1 = strcat(DRIVE, Load_Cell_Paths_Perf2(i));
    SaveFile = sprintf('Perf 2-Wave %d-%d Degrees-Run %d',Wave_Gauge_Limits_Perf2(i,3),...
        Wave_Gauge_Limits_Perf2(i,4), Wave_Gauge_Limits_Perf2(i,5));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] = ...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Perf2(i,1),...
        Wave_Gauge_Limits_Perf2(i,2), SaveFile, SavePath, 2);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(3).T), 0.4375, Wave(3).H);      %Determine Wave Length For Individual Waves
        % Build Nonscalar Structure
        Perf2(n).WaveNum = Wave_Gauge_Limits_Perf2(i,3);
        Perf2(n).Rotation = Wave_Gauge_Limits_Perf2(i,4);
        Perf2(n).RunNum = Wave_Gauge_Limits_Perf2(i,5);
        Perf2(n).GaugeNum = k;
        Perf2(n).WaveTime = WavesTime(:);
        Perf2(n).WaveLimits = Wave(k).Waves;
        Perf2(n).WaveEta = Waves(:,k);
        Perf2(n).WaveHeight = AvgH(k);
        Perf2(n).WaveHm0 = Hm0(k);
        Perf2(n).WavePeriod = AvgT(k);
        Perf2(n).WaveLength = L;
        Perf2(n).Ursells = AvgH(k)/100*L^2/(0.4375^3);
        Perf2_Loc(n,:) = {Perf2(n).WaveNum, Perf2(n).Rotation,...
            Perf2(n).RunNum, Perf2(n).GaugeNum, n};
        n = n+1;
    end

    % Determine max force
    [ForceSeries, AveragedForces, MaxForce] = ...
    Load_Cell_Force(Path1, Load_Cell_Limits_Perf2(i,1), ...
    Load_Cell_Limits_Perf2(i,2), Load_Cell_Limits_Perf2(i,3), ...
    Load_Cell_Limits_Perf2(i,4), L2, Wave(3).H, SaveFile, SavePath);


    % Create structure to hold force time series forces averaged and
    % max force for a given wave and run.
    Perf2_Force(i).WaveNum = Load_Cell_Limits_Perf2(i,5);
    Perf2_Force(i).RunNum = Load_Cell_Limits_Perf2(i,7);
    Perf2_Force(i).Rotation = Load_Cell_Limits_Perf2(i,6);
    Perf2_Force(i).ForceSeries = ForceSeries;
    Perf2_Force(i).AveragedForces = AveragedForces;
    Perf2_Force(i).MaxForce = MaxForce;

    clear WavesTime
    clear Waves
    clear ForceSeries
    clear AveragedForces
    clear MaxForce
end

%% Perf3 Wave Gauge Data Load In
%Initialize locations table
sz = [1,5];
Vartypes = ["double", "double", "double", "double", "double"];
Varnames = ["Wave #", "Run #", "Rotation", "Gauge #", "Location"];
Perf3_Loc = table('Size', sz, 'VariableTypes',Vartypes,'VariableNames',Varnames);

%Initialize location
n = 1;

%Loop through paths and assign data to a nonscalar stucture
for i = 1:length(Wave_Gauge_Paths_Perf3)
    %Initialize Variables
    Path = strcat(DRIVE,Wave_Gauge_Paths_Perf3(i,1));
    Path1 = strcat(DRIVE, Load_Cell_Paths_Perf3(i));
    SaveFile = sprintf('Perf 3-Wave %d-%d Degrees-Run %d',Wave_Gauge_Limits_Perf3(i,3),...
        Wave_Gauge_Limits_Perf3(i,4), Wave_Gauge_Limits_Perf3(i,5));
    SavePath = strcat(SCALED_TESTS_ROOT_DIR, "Figures\", SaveFile);
    
    %Load in wave data using Three_Gauge_Wave_Data function
    [WavesTime(:,1), Waves(:,1:3), AvgH(:), AvgT(:), Wave, Hm0] =...
        Three_Gauge_Wave_Data(Path, Wave_Gauge_Limits_Perf3(i,1),...
        Wave_Gauge_Limits_Perf3(i,2), SaveFile, SavePath);
    for k = 1:3 % Loop Through Each Gauge
        [~,L] = WaveSolver(AvgT(k), 0.4375, AvgH(k)); %Determine Wave Length
        [~,L2] = WaveSolver(seconds(Wave(3).T), 0.4375, Wave(3).H);      %Determine Wave Length For Individual Waves
        % Build Nonscalar Structure
        Perf3(n).WaveNum = Wave_Gauge_Limits_Perf3(i,3);
        Perf3(n).Rotation = Wave_Gauge_Limits_Perf3(i,4);
        Perf3(n).RunNum = Wave_Gauge_Limits_Perf3(i,5);
        Perf3(n).GaugeNum = k;
        Perf3(n).WaveTime = WavesTime(:);
        Perf3(n).WaveLimits = Wave(k).Waves;
        Perf3(n).WaveEta = Waves(:,k);
        Perf3(n).WaveHeight = AvgH(k);
        Perf3(n).WaveHm0 = Hm0(k);
        Perf3(n).WavePeriod = AvgT(k);
        Perf3(n).WaveLength = L;
        Perf3(n).Ursells = AvgH(k)/100*L^2/(0.4375^3);
        Perf3_Loc(n,:) = {Perf3(n).WaveNum, Perf3(n).Rotation,...
            Perf3(n).RunNum, Perf3(n).GaugeNum, n};
        n = n+1;
    end

    % Determine max force
    [ForceSeries, AveragedForces, MaxForce] = ...
    Load_Cell_Force(Path1, Load_Cell_Limits_Perf3(i,1), ...
    Load_Cell_Limits_Perf3(i,2), Load_Cell_Limits_Perf3(i,3), ...
    Load_Cell_Limits_Perf3(i,4), L2, Wave(3).H, SaveFile, SavePath);


    % Create structure to hold force time series forces averaged and
    % max force for a given wave and run.
    Perf3_Force(i).WaveNum = Load_Cell_Limits_Perf3(i,5);
    Perf3_Force(i).RunNum = Load_Cell_Limits_Perf3(i,7);
    Perf3_Force(i).Rotation = Load_Cell_Limits_Perf3(i,6);
    Perf3_Force(i).ForceSeries = ForceSeries;
    Perf3_Force(i).AveragedForces = AveragedForces;
    Perf3_Force(i).MaxForce = MaxForce;

    clear WavesTime
    clear Waves
    clear ForceSeries
    clear AveragedForces
    clear MaxForce
end

%% Force Averaging for each perforation and wave

%Ref Monopile
TempWave1 = [];
TempWave2 = [];
TempWave3 = [];

for l = 1:numel(Ref_Force)
    if Ref_Force(l).WaveNum == 1
        n = length(TempWave1)+1;
        TempWave1(n) = Ref_Force(l).MaxForce;
    elseif Ref_Force(l).WaveNum == 2
        if Ref_Force(l).MaxForce < 14
            n = length(TempWave2)+1;
            TempWave2(n) = Ref_Force(l).MaxForce;
        end
    else
        n = length(TempWave3)+1;
        TempWave3(n) = Ref_Force(l).MaxForce;
    end
end
Ref_Avg_Force = [mean(TempWave1), mean(TempWave2), mean(TempWave3)];


%Perf 1 Monopile
TempArray = [Perf1_Force.WaveNum; Perf1_Force.Rotation; Perf1_Force.MaxForce];
TempArray = transpose(TempArray);
TempArray1 = TempArray(TempArray(:,1)==1,2:3);
TempArray2 = TempArray(TempArray(:,1)==2,2:3);
TempArray3 = TempArray(TempArray(:,1)==3,2:3);
TempWave1 = [mean(TempArray1(TempArray1(:,1)==0,2)), mean(TempArray1(TempArray1(:,1)==15,2)),...
    mean(TempArray1(TempArray1(:,1)==30,2)), mean(TempArray1(TempArray1(:,1)==45,2))];
TempWave2 = [mean(TempArray2(TempArray2(:,1)==0,2)), mean(TempArray2(TempArray2(:,1)==15,2)),...
    mean(TempArray2(TempArray2(:,1)==30,2)), mean(TempArray2(TempArray2(:,1)==45,2))];
TempWave3 = [mean(TempArray3(TempArray3(:,1)==0,2)), mean(TempArray3(TempArray3(:,1)==15,2)),...
    mean(TempArray3(TempArray3(:,1)==30,2)), mean(TempArray3(TempArray3(:,1)==45,2))];
Perf1_Avg_Forces = [TempWave1; TempWave2; TempWave3];


%Perf 2 Monopile
TempArray = [Perf2_Force.WaveNum; Perf2_Force.Rotation; Perf2_Force.MaxForce];
TempArray = transpose(TempArray);
TempArray1 = TempArray(TempArray(:,1)==1,2:3);
TempArray2 = TempArray(TempArray(:,1)==2,2:3);
TempArray3 = TempArray(TempArray(:,1)==3,2:3);
TempWave1 = [mean(TempArray1(TempArray1(:,1)==0,2)), mean(TempArray1(TempArray1(:,1)==15,2)),...
    mean(TempArray1(TempArray1(:,1)==30,2)), mean(TempArray1(TempArray1(:,1)==45,2))];
TempWave2 = [mean(TempArray2(TempArray2(:,1)==0,2)), mean(TempArray2(TempArray2(:,1)==15,2)),...
    mean(TempArray2(TempArray2(:,1)==30,2)), mean(TempArray2(TempArray2(:,1)==45,2))];
TempWave3 = [mean(TempArray3(TempArray3(:,1)==0,2)), mean(TempArray3(TempArray3(:,1)==15,2)),...
    mean(TempArray3(TempArray3(:,1)==30,2)), mean(TempArray3(TempArray3(:,1)==45,2))];
Perf2_Avg_Forces = [TempWave1; TempWave2; TempWave3];


%Perf 3 Monopile
TempArray = [Perf3_Force.WaveNum; Perf3_Force.Rotation; Perf3_Force.MaxForce];
TempArray = transpose(TempArray);
TempArray1 = TempArray(TempArray(:,1)==1,2:3);
TempArray2 = TempArray(TempArray(:,1)==2,2:3);
TempArray3 = TempArray(TempArray(:,1)==3,2:3);
TempWave1 = [mean(TempArray1(TempArray1(:,1)==0,2)), mean(TempArray1(TempArray1(:,1)==15,2)),...
    mean(TempArray1(TempArray1(:,1)==30,2)), mean(TempArray1(TempArray1(:,1)==45,2))];
TempWave2 = [mean(TempArray2(TempArray2(:,1)==0,2)), mean(TempArray2(TempArray2(:,1)==15,2)),...
    mean(TempArray2(TempArray2(:,1)==30,2)), mean(TempArray2(TempArray2(:,1)==45,2))];
TempWave3 = [mean(TempArray3(TempArray3(:,1)==0,2)), mean(TempArray3(TempArray3(:,1)==15,2)),...
    mean(TempArray3(TempArray3(:,1)==30,2)), mean(TempArray3(TempArray3(:,1)==45,2))];
Perf3_Avg_Forces = [TempWave1; TempWave2; TempWave3];

%% Determine comparison values 
% Perforation 1
Perf1_Rel_Force = [Perf1_Avg_Forces(1,:)./Ref_Avg_Force(1);...
    Perf1_Avg_Forces(2,:)./Ref_Avg_Force(2);...
    Perf1_Avg_Forces(3,:)./Ref_Avg_Force(3)];

% Perforation 2
Perf2_Rel_Force = [Perf2_Avg_Forces(1,:)./Ref_Avg_Force(1);...
    Perf2_Avg_Forces(2,:)./Ref_Avg_Force(2);...
    Perf2_Avg_Forces(3,:)./Ref_Avg_Force(3)];

% Perforation 3
Perf3_Rel_Force = [Perf3_Avg_Forces(1,:)./Ref_Avg_Force(1);...
    Perf3_Avg_Forces(2,:)./Ref_Avg_Force(2);...
    Perf3_Avg_Forces(3,:)./Ref_Avg_Force(3)];

%% Lines of Best Fit for scaled CFD data
x = [KC(1), KC(1), KC(1), KC(1),KC(3), KC(3), KC(3), KC(3), KC(2), KC(2), KC(2), KC(2)];
s2(1,:) = tanhfit(x,[Perf1_Rel_Force(1,:), Perf1_Rel_Force(3,:), Perf1_Rel_Force(2,:)]);
s2(2,:) = tanhfit(x,[Perf2_Rel_Force(1,:), Perf2_Rel_Force(3,:), Perf2_Rel_Force(2,:)]);
s2(3,:) = tanhfit(x,[Perf3_Rel_Force(1,:), Perf3_Rel_Force(3,:), Perf3_Rel_Force(2,:)]);

%% Root Mean Square Error (RMSE) in relation to other research
s = [0.9857, 0.9238; 0.9092, 0.6742; 0.6764, 0.6183];
x1 = linspace(0, 10);
bestfit=@(b1,b2,x) b1.*tanh(b2.*x);
RMSE_Ploeg(1) = sqrt(sum(bestfit(1,1.49,x1) - bestfit(s2(1,1),s2(1,2),x1)).^2)/length(x1);
RMSE_Ploeg(2) = sqrt(sum(bestfit(1,1.49,x1) - bestfit(s2(2,1),s2(2,2),x1)).^2)/length(x1);
RMSE_Ploeg(3) = sqrt(sum(bestfit(1,1.49,x1) - bestfit(s2(3,1),s2(3,2),x1)).^2)/length(x1);

RMSE_CFD(1) = sqrt(sum(bestfit(s(1,1),s(1,2),x1) - bestfit(s2(1,1),s2(1,2),x1)).^2)/length(x1);
RMSE_CFD(2) = sqrt(sum(bestfit(s(2,1),s(2,2),x1) - bestfit(s2(2,1),s2(2,2),x1)).^2)/length(x1);
RMSE_CFD(3) = sqrt(sum(bestfit(s(3,1),s(3,2),x1) - bestfit(s2(3,1),s2(3,2),x1)).^2)/length(x1);

RMSE_Ander(1) = sqrt(sum((bestfit(.94,0.84,x1) - bestfit(s2(2,1),s2(2,2),x1)).^2)/length(x1));


%% R-Squared Calculation
R_SQ_Ploeg = [R_SQ(bestfit(s2(1,1),s2(1,2),x1), bestfit(1,1.49,x1));...
    R_SQ(bestfit(s2(2,1),s2(2,2),x1), bestfit(0.93,0.75,x1));...
    R_SQ(bestfit(s2(3,1),s2(3,2),x1), bestfit(0.75,0.52,x1))];
R_SQ_CFD = [R_SQ(bestfit(s2(1,1),s2(1,2),x1), bestfit(s(1,1),s(1,2),x1));...
    R_SQ(bestfit(s2(2,1),s2(2,2),x1), bestfit(s(2,1),s(2,2),x1));...
    R_SQ(bestfit(s2(3,1),s2(3,2),x1), bestfit(s(3,1),s(3,2),x1))];
R_SQ_Andersen = R_SQ(bestfit(s2(2,1),s2(2,2),x1), bestfit(0.94,0.84,x1));

%% Plotting Time


% Perforation 1
figure(1)
plot(KC,Perf1_Rel_Force(:,1),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf1_Rel_Force(:,2),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf1_Rel_Force(:,3),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Perf1_Rel_Force(:,4),"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 1*tanh(1.49*x), [1 6], 'g');
fplot(@(x) 0.99*tanh(0.92*x), [1 6], 'r');
fplot(@(x) s2(1,1)*tanh(s2(1,2)*x), [1 6], 'm');
ylim([0.3 1.1])
title("Perforation #1 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg (2021)','CFD Line of Best Fit', "WT Line of Best Fit")
functionText1 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(1));
functionText2 = sprintf('CFD = %.4f', R_SQ_CFD(1));
text(1.25, .45, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(1.5, .41, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');

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
fplot(@(x) s2(2,1)*tanh(s2(2,2)*x), [1 6], 'm');
ylim([0.3 1.1])
title("Perforation #2 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', ...
    'Anderson, et al. (2020)', 'Ploeg (2021)','CFD Line of Best Fit', 'WT Line of Best Fit')
functionText1 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(2));
functionText2 = sprintf('Andersen, et al. (2020) = %.4f', R_SQ_Andersen(1));
functionText3 = sprintf('CFD = %.4f', R_SQ_CFD(2));
text(1.25, .45, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(1.5, .41, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(1.5, .37, functionText3, 'FontSize', 10, 'HorizontalAlignment', 'left');


% Perforation 3
figure(3)
plot(KC,Perf3_Rel_Force(:,1),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf3_Rel_Force(:,2),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf3_Rel_Force(:,3),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
plot(KC,Perf3_Rel_Force(:,4),"Color",[1 1 0],"Marker","square","LineStyle","none")
fplot(@(x) 0.75*tanh(0.52*x), [1 6], 'g');
fplot(@(x) 0.68*tanh(0.62*x), [1 6], 'r');
fplot(@(x) s2(3,1)*tanh(s2(3,2)*x), [1 6], 'm');
ylim([0.3 1.1])
title("Perforation #3 Relative Force vs. KC number")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
legend('0 Degrees', '15 Degrees', '30 Degrees', '45 Degrees', 'Ploeg (2021)','CFD Line of Best Fit', 'WT Line of Best Fit')
functionText1 = sprintf('R^2: Ploeg (2021) = %.4f', R_SQ_Ploeg(3));
functionText2 = sprintf('CFD = %.4f', R_SQ_CFD(3));
text(2, .45, functionText1, 'FontSize', 10, 'HorizontalAlignment', 'left');
text(2.25, .41, functionText2, 'FontSize', 10, 'HorizontalAlignment', 'left');





%% Create table output of all max force values
rownames = {'Reference','Perforation 1', 'Perforation 2', 'Perforation 3'};
columnnames = {'Wave 1', 'Wave 2', 'Wave3'};
Wave1Peaks =[Ref_Avg_Force(1); mean(Perf1_Avg_Forces(1,:)); mean(Perf2_Avg_Forces(1,:)); mean(Perf3_Avg_Forces(1,:))];
Wave2Peaks =[Ref_Avg_Force(2); mean(Perf1_Avg_Forces(2,:)); mean(Perf2_Avg_Forces(2,:)); mean(Perf3_Avg_Forces(2,:))];
Wave3Peaks =[Ref_Avg_Force(3); mean(Perf1_Avg_Forces(3,:)); mean(Perf2_Avg_Forces(3,:)); mean(Perf3_Avg_Forces(3,:))];
T1 = table(Wave1Peaks, Wave2Peaks, Wave3Peaks, 'RowNames',rownames, 'VariableNames',columnnames)

%% Create table output of average wave parameters
RefWaveParam = [extractfield(Ref,'WaveHeight'); extractfield(Ref,'WavePeriod'); extractfield(Ref,'WaveLength')];
Perf1WaveParam = [extractfield(Perf1,'WaveHeight'); extractfield(Perf1,'WavePeriod'); extractfield(Perf1,'WaveLength')];
Perf2WaveParam = [extractfield(Perf2,'WaveHeight'); extractfield(Perf2,'WavePeriod'); extractfield(Perf2,'WaveLength')];
Perf3WaveParam = [extractfield(Perf3,'WaveHeight'); extractfield(Perf3,'WavePeriod'); extractfield(Perf3,'WaveLength')];

rownames = {'Height','Period', 'Length'};
columnnames = {'Wave 1', 'Wave 2', 'Wave3'};
WaveGauge3 = [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108];
Wave1 =[mean([RefWaveParam(1,WaveGauge3(1:3)), Perf1WaveParam(1,WaveGauge3(1:12)), Perf2WaveParam(1,WaveGauge3(1:12)), Perf3WaveParam(1,WaveGauge3(1:11))]);...
    mean([RefWaveParam(2,WaveGauge3(1:3)), Perf1WaveParam(2,WaveGauge3(1:12)), Perf2WaveParam(2,WaveGauge3(1:12)), Perf3WaveParam(2,WaveGauge3(1:11))]); ...
    mean([RefWaveParam(3,WaveGauge3(1:3)), Perf1WaveParam(3,WaveGauge3(1:12)), Perf2WaveParam(3,WaveGauge3(1:12)), Perf3WaveParam(3,WaveGauge3(1:11))])];
Wave2 =[mean([RefWaveParam(1,WaveGauge3(4:6)), Perf1WaveParam(1,WaveGauge3(13:24)), Perf2WaveParam(1,WaveGauge3(13:24)), Perf3WaveParam(1,WaveGauge3(12:22))]);...
    mean([RefWaveParam(2,WaveGauge3(4:6)), Perf1WaveParam(2,WaveGauge3(13:24)), Perf2WaveParam(2,WaveGauge3(13:24)), Perf3WaveParam(2,WaveGauge3(12:22))]);...
    mean([RefWaveParam(3,WaveGauge3(4:6)), Perf1WaveParam(3,WaveGauge3(13:24)), Perf2WaveParam(3,WaveGauge3(13:24)), Perf3WaveParam(3,WaveGauge3(12:22))])];
Wave3 =[mean([RefWaveParam(1,WaveGauge3(7:9)), Perf1WaveParam(1,WaveGauge3(25:36)), Perf2WaveParam(1,WaveGauge3(25:36)), Perf3WaveParam(1,WaveGauge3(23:33))]);...
    mean([RefWaveParam(2,WaveGauge3(7:9)), Perf1WaveParam(2,WaveGauge3(25:36)), Perf2WaveParam(2,WaveGauge3(25:36)), Perf3WaveParam(2,WaveGauge3(23:33))]);...
    mean([RefWaveParam(3,WaveGauge3(7:9)), Perf1WaveParam(3,WaveGauge3(25:36)), Perf2WaveParam(3,WaveGauge3(25:36)), Perf3WaveParam(3,WaveGauge3(23:33))])];
T2 = table(Wave1, Wave2, Wave3, 'RowNames',rownames, 'VariableNames',columnnames)

%% All Perforations Combined
Perf1_Avg_Rel_Force = [mean(Perf1_Rel_Force(1,:)), mean(Perf1_Rel_Force(2,:)), mean(Perf1_Rel_Force(3,:))];
Perf2_Avg_Rel_Force = [mean(Perf2_Rel_Force(1,:)), mean(Perf2_Rel_Force(2,:)), mean(Perf2_Rel_Force(3,:))];
Perf3_Avg_Rel_Force = [mean(Perf3_Rel_Force(1,:)), mean(Perf3_Rel_Force(2,:)), mean(Perf3_Rel_Force(3,:))];

figure(4)
plot(KC,Perf1_Avg_Rel_Force(:),"Color",[1 0 0],"Marker","+","LineStyle","none")
hold on
plot(KC,Perf2_Avg_Rel_Force(:),"Color",[0 1 0],"Marker","o","LineStyle","none")
plot(KC,Perf3_Avg_Rel_Force(:),"Color",[0 0 1],"Marker","diamond","LineStyle","none")
fplot(@(x) s2(1,1)*tanh(s2(1,2)*x), [1 6], 'r');
fplot(@(x) s2(2,1)*tanh(s2(2,2)*x), [1 6], 'g');
fplot(@(x) s2(3,1)*tanh(s2(3,2)*x), [1 6], 'b');
ylim([0.3 1.1])
title("Relative Force vs. KC number - Wave Tank Results")
xlabel("KC")
ylabel("Perforation Force/Reference Force")
legend('Perforation 1', 'Perforation 2', 'Perforation 3', 'Perforation 1 Line of Best Fit', 'Perforation 2 Line of Best Fit','Perforation 3 Line of Best Fit')