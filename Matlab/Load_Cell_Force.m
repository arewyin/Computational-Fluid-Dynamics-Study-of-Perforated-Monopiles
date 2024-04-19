function [Force, Force_Maxs, Force_Max] = Load_Cell_Force(data_file_path,...
    zero_start, zero_end, data_start, data_end, L, H, SaveFile, SavePath)
%Inputs strain from a given file and outputs a series of forces using the
%calibration coefficient determined by Load_Cell_Calb2.m
%Input Variables:
%       data_file_path = string location of file to load
%       zero_start = Double of start of zeroed data
%       zero_end = Double of end of zeroed data
%       data_start = start of windowed data
%       data_end = end of windowed data
%       L = Array of Coordinating Wave lengths
%       H = Array of Coordinating Wave Heights


%% Import Data
Strain = table2array(importfile2(data_file_path, [1, Inf]));

%% Correct Input Variables
if data_end == 0
    data_end = length(Strain);
end

if zero_end == 0
    zero_end = length(Strain);
end


%% Determine Zero of Signal
Average_zero = mean(Strain(zero_start:zero_end));

%% Determine "Waves" in zero strain data
ZeroedStrain = Strain - Average_zero;
[IndexWaves] = StrainZeroCrossings(ZeroedStrain(data_start:data_end));

%% Limits
if length(L) > length(IndexWaves)-1
    Limit = length(IndexWaves)-1;
else
    Limit = length(L);
end

%% Determine Force
Force = [];
for i = 1:Limit
    [ZLoca(i).Centroids ZLoca(i).Force] = MFCentroid(H(i), L(i));
    z = @(x,y) 4.16681*10^-5*exp(1.48429*ZLoca(x).Centroids(y));
    StartIndex = IndexWaves(i) + data_start-1;
    EndIndex = IndexWaves(i+1) + data_start-2;
    WaveStrain = ZeroedStrain(StartIndex:EndIndex);
    [MaxStrain, MaxStrainIndex] = max(WaveStrain); 
    [MinStrain, MinStrainIndex] = min(WaveStrain);
    if MaxStrainIndex > MinStrainIndex %Checks if Crest is after Trough
        Condition = 1;

        %Deterimining all values above 0
        surzero = find(ZLoca(i).Force > 0);

        %% Check for Downcrossing
        n=0;
        %Determine limit of for loop
        if surzero(length(surzero)) == length(ZLoca(i).Force)
            limit = length(surzero)-1;
        else
            limit = length(surzero);
        end

        %For loop to check for down crossings
        for l = 1:limit
            if ZLoca(i).Force(surzero(l)+1) <= 0 %Check if next value is below zero
            n=n+1;
            Waves(n) = surzero(l) + 1;
            end
        end

        Mid1 = MinStrainIndex;
        Mid2 = MaxStrainIndex;
        [CrestHeight(i), Point2] = max(ZLoca(i).Force(Waves(1):Waves(2)));
        [TroughHeight(i), Point1] = min(ZLoca(i).Force(Waves(1):Waves(2)));

    elseif MaxStrainIndex < MinStrainIndex %Checks if Crest is before Trough
        Condition = 2;

        %Deterimining all values below 0
        subzero = find(ZLoca(i).Force < 0);

        %Check for Upcrossing
        n=0;
        %Determine limit of for loop
        if length(subzero) == length(ZLoca(i).Force)
            limit = length(subzero)-1;
        else
            limit = length(subzero);
        end 

        for l = 1:limit
            if ZLoca(i).Force(subzero(l)+1) >= 0 %Check if next value is above zero
                n=n+1;
                Waves(n) = subzero(l) + 1;
            end
        end

        Mid1 = MaxStrainIndex;
        Mid2 = MinStrainIndex;
        [CrestHeight(i), Point1] = max(ZLoca(i).Force(Waves(1):Waves(2)));
        [TroughHeight(i), Point2] = min(ZLoca(i).Force(Waves(1):Waves(2)));
    end

    % Correcting Points
    Point1 = Point1 + Waves(1);
    Point2 = Point2 + Waves(1);

    % Add Zeros to WaveStrain if it doesn't end already
    if WaveStrain(1) ~= 0
        WaveStrain = [0; WaveStrain];
    end
    if WaveStrain(end) ~= 0 
        WaveStrain = [WaveStrain; 0];
    end

    %Expand Strain Data to match length of force data aligning peaks,
    % troughs and zeros
    xq1 = linspace(1,length(WaveStrain(1:Mid1+1)),abs(Waves(1) - Point1));
    x1 = 1:1:length(WaveStrain(1:Mid1+1));
    yq1 = interp1(x1, WaveStrain(1:Mid1+1), xq1);

    xq2 = linspace(1,length(WaveStrain(Mid1+1:Mid2+1)),abs(Point1 - Point2));
    x2 = 1:1:length(WaveStrain(Mid1+1:Mid2+1));
    yq2 = interp1(x2, WaveStrain(Mid1+1:Mid2+1), xq2);

    xq3 = linspace(1,length(WaveStrain(Mid2+1:end)),abs(Point2-Waves(2)-1));
    x3 = 1:1:length(WaveStrain(Mid2+1:end));
    yq3 = interp1(x3, WaveStrain(Mid2+1:end), xq3);

    yq(i,:) = [yq1, yq2, yq3];

    for j = 1:length(yq(i,:))
        TempForce(j) = yq(i,j)*z(i,j+Waves(1));
    end
    Force = [Force,TempForce];
end

%% Max Wave Force
n = Limit*2;
Force_Maxs(:) = maxk(abs(Force),n);
Force_Max = mean(Force_Maxs);

% SavePath = strcat(SavePath,"_LoadCell");
% SaveFile = strcat(SaveFile,"-Force");
% fig = figure();
% plot(Force)
% title(SaveFile + '-Force')
% ylabel('Force (N)')
% xlabel('Sample Number')
% saveas(fig,SavePath,'png');
% close()

end