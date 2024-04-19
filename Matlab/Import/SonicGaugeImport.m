function Test = SonicGaugeImport(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  TEST1000 = IMPORTFILE3(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  TEST1000 = IMPORTFILE3(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  Test = SonicGaugeImport("F:\Thesis\Wave Tank Tests\Wave Gauge Tests\Sonic Gauges\20231026-Test1000.csv", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 26-Oct-2023 16:00:18

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 36);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["F003", "S32", "B0", "T101107", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36"];
opts.VariableTypes = ["double", "double", "double", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "T101107", "InputFormat", "uuuu-MM-dd'T'HH:mm:ss");
opts = setvaropts(opts, ["F003", "S32", "B0"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["F003", "S32", "B0"], "ThousandsSeparator", ",");

% Import the data
Test = readtable(filename, opts);

end