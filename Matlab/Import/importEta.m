function Perforation110Eta = importEta(filename, dataLines)
%IMPORTFILE3 Import data from a text file
%  PERFORATION110ETA = IMPORTFILE3(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  PERFORATION110ETA = IMPORTFILE3(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Perforation1_1_0-Eta = importEta("F:\Thesis\Scaled Simulations\Wave Tank Parameters\Perforation_1_1_0_Eta.txt", [10, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 28-Feb-2024 16:06:52

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [10, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["time", "hyel_1", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"];
opts.SelectedVariableNames = ["time", "hyel_1"];
opts.VariableTypes = ["double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "time", "TrimNonNumeric", true);
opts = setvaropts(opts, "time", "ThousandsSeparator", ",");

% Import the data
Perforation110Eta = readtable(filename, opts);

end