function WaveGaugePathsRef = importRefPaths(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  WAVEGAUGEPATHSREF = IMPORTFILE3(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  WAVEGAUGEPATHSREF = IMPORTFILE3(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  WaveGaugePathsRef = importRefPaths("D:\Thesis\Wave Tank Tests\Scaled Tests\Wave_Gauge_Paths_Ref.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 22-Jan-2024 16:15:16

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Path", "Start", "End", "Wave", "Run"];
opts.VariableTypes = ["string", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Path", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Path", "EmptyFieldRule", "auto");

% Import the data
WaveGaugePathsRef = readtable(filename, opts);

end