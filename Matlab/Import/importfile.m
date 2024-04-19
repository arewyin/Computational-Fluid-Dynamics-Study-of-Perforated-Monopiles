function Reference13Data = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  REFERENCE13DATA = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  REFERENCE13DATA = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Reference13Data = importfile("F:\School\FIT\Thesis\Scaled Simulations\Reference_13-Data.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 23-Aug-2023 15:33:49

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 147);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Component1CombinedFluidForceMagnitude", "Component1CombinedFluidTorqueMagnitude", "Component1PressureForce", "Component1ShearForce", "Component1Volume", "Component1XCombinedFluidForce", "Component1XCombinedFluidTorqueAboutReferencePoint", "Component1XPressureForce", "Component1XShearForce", "Component1YCombinedFluidForce", "Component1YCombinedFluidTorqueAboutReferencePoint", "Component1YPressureForce", "Component1YShearForce", "Component1ZCombinedFluidForce", "Component1ZCombinedFluidTorqueAboutReferencePoint", "Component1ZPressureForce", "Component1ZShearForce", "DiagnosticsParallelEfficiency", "DiagnosticsAccumulatedFluidMassSource", "DiagnosticsConvectiveVolumeError", "DiagnosticsConvectiveVolumeErrorLost", "DiagnosticsElapsedClockTime", "DiagnosticsElapsedClockTimePerTimeStep", "DiagnosticsFluid1VolumeNetInflux", "DiagnosticsFreeSurfaceTimestepLimit", "DiagnosticsInterblockBoundaryVolumeError", "DiagnosticsInterblockBoundaryVolumeErrorLost", "DiagnosticsMaximumPressureResidual", "DiagnosticsPressureConvergenceCriterion", "DiagnosticsPressureIterationCount", "DiagnosticsPressureRelaxationFactor", "DiagnosticsTimestepSize", "DiagnosticsTimestepStabilityLimit", "DiagnosticsViscousTimestepLimit", "DiagnosticsVolumetricConvectiveTimestepLimit", "DiagnosticsXdirectionConvectiveTimestepLimit", "DiagnosticsYdirectionConvectiveTimestepLimit", "DiagnosticsZdirectionConvectiveTimestepLimit", "F3DGlobal_0", "F3DGlobal_1", "F3DGlobal_2", "F3DGlobal_3", "F3DGlobal_4", "F3DGlobal_5", "F3DGlobal_6", "F3DGlobal_Magnitude", "F3DProcess_0", "F3DProcess_1", "F3DProcess_2", "F3DProcess_Magnitude", "F3D_COMPONENT_FACETS", "F3D_COMPONENT_TYPES", "F3D_FIXED_DATA_0", "F3D_FIXED_DATA_1", "F3D_HISTORY_DATA_0", "F3D_HISTORY_DATA_1", "F3D_HISTORY_DATA_2", "F3D_MESH_ARRAY_TYPES", "F3D_MESH_DATA_0", "F3D_MESH_DATA_1", "F3D_PARTICLE_DATA_0", "F3D_PARTICLE_DATA_1", "FillFraction", "FluidCenterOfMassXcoordinate", "FluidCenterOfMassYcoordinate", "FluidCenterOfMassZcoordinate", "FluidSurfaceArea", "HistoryProbe1DynamicViscosity", "HistoryProbe1FluidVelocityMagnitude", "HistoryProbe1FluidXvelocity", "HistoryProbe1FluidYvelocity", "HistoryProbe1FluidZvelocity", "HistoryProbe1FractionOfFluid", "HistoryProbe1Pressure", "HistoryProbe1ProbeVelocityMagnitude", "HistoryProbe1ProbeXvelocity", "HistoryProbe1ProbeYvelocity", "HistoryProbe1ProbeZvelocity", "HistoryProbe1StrainRateMagnitude", "HistoryProbe1TurbulenceIntensity", "HistoryProbe1TurbulentDissipation", "HistoryProbe1TurbulentEnergy", "HistoryProbe1Xcoordinate", "HistoryProbe1Ycoordinate", "HistoryProbe1Zcoordinate", "MassaveragedFluidMeanKineticEnergy", "MassaveragedTurbulentBuoyantProduction", "MassaveragedTurbulentDissipation", "MassaveragedTurbulentKineticEnergy", "MassaveragedTurbulentShearProduction", "MaximumFluid1Velocity", "MaximumFluidSurfaceVelocity", "MeshBlock1", "YmaxFluid1VolumeFlowRate", "MeshBlock2", "YminFluid1Elevation", "MeshBlock3", "YminFluid1VolumeFlowRate", "MeshBlock4", "YminSpecifiedTurbulentDissipation", "MeshBlock5", "YminSpecifiedTurbulentEnergy", "MeshBlock6", "YminSpecifiedXvelocityOfCurrent", "MeshBlock7", "YminSpecifiedYvelocityOfCurrent", "MeshBlock8", "YminSpecifiedZvelocityOfCurrent", "MeshBlock9", "ZmaxFluid1VolumeFlowRate", "MeshBlock10", "ZmaxSpecifiedFluidFraction", "MeshBlock11", "ZmaxSpecifiedPressure", "MeshBlock12", "ZmaxSpecifiedTurbulentDissipation", "MeshBlock13", "ZmaxSpecifiedTurbulentEnergy", "MeshBlock14", "ZmaxSpecifiedXvelocity", "MeshBlock15", "ZmaxSpecifiedYvelocity", "MeshBlock16", "XmaxFluid1VolumeFlowRate", "MeshBlock17", "XminFluid1VolumeFlowRate", "MeshBlock18", "YmaxFluid1VolumeFlowRate1", "MeshBlock19", "YminFluid1VolumeFlowRate1", "MeshBlock20", "ZmaxFluid1VolumeFlowRate1", "MeshBlock21", "ZmaxSpecifiedFluidFraction1", "MeshBlock22", "ZmaxSpecifiedPressure1", "MeshBlock23", "ZmaxSpecifiedTurbulentDissipation1", "MeshBlock24", "ZmaxSpecifiedTurbulentEnergy1", "MeshBlock25", "ZmaxSpecifiedXvelocity1", "MeshBlock26", "ZmaxSpecifiedYvelocity1", "NumberOfMeshBlocks", "ProbeParticles", "Time", "TotalFluidMass", "TotalNumberOfParticles", "VoidTotalVolume", "VoidVolumeCenterX", "VoidVolumeCenterY", "VoidVolumeCenterZ", "VolumeOfFluid1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "categorical", "categorical", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["F3D_FIXED_DATA_0", "F3D_FIXED_DATA_1", "F3D_HISTORY_DATA_0", "F3D_HISTORY_DATA_1", "F3D_MESH_DATA_0", "F3D_MESH_DATA_1", "F3D_PARTICLE_DATA_0", "F3D_PARTICLE_DATA_1", "YminFluid1VolumeFlowRate1", "MeshBlock20", "ZmaxFluid1VolumeFlowRate1", "MeshBlock21", "ZmaxSpecifiedFluidFraction1", "MeshBlock22", "ZmaxSpecifiedPressure1", "MeshBlock23", "ZmaxSpecifiedTurbulentDissipation1", "MeshBlock24", "ZmaxSpecifiedTurbulentEnergy1", "MeshBlock25", "ZmaxSpecifiedXvelocity1", "MeshBlock26", "ZmaxSpecifiedYvelocity1", "NumberOfMeshBlocks", "ProbeParticles", "Time", "TotalFluidMass", "TotalNumberOfParticles", "VoidTotalVolume", "VoidVolumeCenterX", "VoidVolumeCenterY", "VoidVolumeCenterZ", "VolumeOfFluid1"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["F3D_FIXED_DATA_0", "F3D_FIXED_DATA_1", "F3D_HISTORY_DATA_0", "F3D_HISTORY_DATA_1", "F3D_HISTORY_DATA_2", "F3D_MESH_ARRAY_TYPES", "F3D_MESH_DATA_0", "F3D_MESH_DATA_1", "F3D_PARTICLE_DATA_0", "F3D_PARTICLE_DATA_1", "YminFluid1VolumeFlowRate1", "MeshBlock20", "ZmaxFluid1VolumeFlowRate1", "MeshBlock21", "ZmaxSpecifiedFluidFraction1", "MeshBlock22", "ZmaxSpecifiedPressure1", "MeshBlock23", "ZmaxSpecifiedTurbulentDissipation1", "MeshBlock24", "ZmaxSpecifiedTurbulentEnergy1", "MeshBlock25", "ZmaxSpecifiedXvelocity1", "MeshBlock26", "ZmaxSpecifiedYvelocity1", "NumberOfMeshBlocks", "ProbeParticles", "Time", "TotalFluidMass", "TotalNumberOfParticles", "VoidTotalVolume", "VoidVolumeCenterX", "VoidVolumeCenterY", "VoidVolumeCenterZ", "VolumeOfFluid1"], "EmptyFieldRule", "auto");

% Import the data
Reference13Data = readtable(filename, opts);

end