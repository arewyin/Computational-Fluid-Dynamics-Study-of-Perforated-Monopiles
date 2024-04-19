% 3_Gauge_Wave_Data
function [Cap_Time, Cap_Data, AvgH, AvgT, Wave, Hm0] = ...
    Three_Gauge_Wave_Data(data_file_path, Wave_Start, Wave_End, SaveFile, SavePath, k)
% Inputs: data file path, Start of waves, End of waves

% Outputs: Timestamps, Wave Guage data, Average Height of Waves, 
% Average Period of Waves, Array of Individual Wave Heights and an Array of
% Individual Wave Periods

%% Import Data
Cap = CapGauge3Import(data_file_path, [6, Inf]);
Start_Time_Cap = DateTimeCapGauges(data_file_path, [1, 1]);

%% Data Sorting
Cap_Data = [Cap.('Value'), Cap.('Peak'), Cap.('Offset')];
Cap_Time = [Cap.('Time')];

if exist('k','var')
    Cap_Data = Cap_Data*k;
end

%% Fix DateTime
str = string(Start_Time_Cap.('Time'));
d = duration(str);
Start_Time_Cap = Start_Time_Cap.('Date') + d;
Start_Time_Cap = datetime(Start_Time_Cap, "Format","MM-dd-yyyy HH:mm:ss.SSS");
str = string(Cap_Time);
d = duration(str);
Cap_Time = Start_Time_Cap + d;
Cap_Time = datetime(Cap_Time, "Format","MM-dd-yyyy HH:mm:ss.SSS");


%% Convert Data from m to cm
Cap_Data_CM = Cap_Data*100;

%% Determine Wave Parameters
[Wave(1).Waves,Wave(1).H, Wave(1).T, AvgH(1), AvgT(1), ~, ~, ~, ~, ~, ~, Hm0(1)] = ...
    ZeroUp(Cap_Data(Wave_Start:Wave_End,1), Cap_Time(Wave_Start:Wave_End));
[Wave(2).Waves,Wave(2).H, Wave(2).T, AvgH(2), AvgT(2), ~, ~, ~, ~, ~, ~, Hm0(2)] = ...
    ZeroUp(Cap_Data(Wave_Start:Wave_End,2), Cap_Time(Wave_Start:Wave_End));
[Wave(3).Waves,Wave(3).H, Wave(3).T, AvgH(3), AvgT(3), ~, ~, ~, ~, ~, ~, Hm0(3)] = ...
    ZeroUp(Cap_Data(Wave_Start:Wave_End,3), Cap_Time(Wave_Start:Wave_End));

% %% Plot Wave Series for all three gauges
% SavePath = strcat(SavePath,"_WaveGauges");
% 
% fig=figure();
% plot(Cap_Time(Wave_Start:Wave_End),Cap_Data(Wave_Start:Wave_End,1))
% hold on 
% plot(Cap_Time(Wave_Start:Wave_End),Cap_Data(Wave_Start:Wave_End,2))
% plot(Cap_Time(Wave_Start:Wave_End),Cap_Data(Wave_Start:Wave_End,3))
% title(SaveFile)
% xlabel('Time (HH:mm:ss)')
% ylabel('Eta (cm)')
% legend('Gauge #1', 'Gauge #2', 'Gauge #3')
% saveas(fig,SavePath,'png');
% close()
