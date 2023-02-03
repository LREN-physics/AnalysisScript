function Params=GetParams()
% Returns a structure containing the data paths and general parameters  
% regarding the analsyis to be conducted.  
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params.HomeDir='yourpath\Data'; % Analysis output folder
Params.DataDir='yourpath\Data'; % Location of the Data folder

Params.NMDir=fullfile(Params.DataDir,'Neuromorphometrics'); % Location the Neuromorphometrics folder
Params.AnalysisFolder='Analysis';                           % Naming of the analysis folder
Params.DataSubDir='';                                       % Data location within the main folder of each dataset.

Params.AgeModelOrder=3;
Params.MDIvsResOrder=3; % Order of the polynomial fit of var(res)vsMDI

Params.ReMLAnal.RefLambda=3;      % Power of the reference noise model for analysis of the Free Energy change
Params.MotionRegPowers=[1 2 3 4]; % Powers of the MDI inserted in design matrix for MotionBias analyses

% Labelling of grey and white matter regions from the neuromorphometrics
% atlas indices
AtlasVals=[0;4;11;23;30;31;32;35;36;37;38;39;40;41;44;45;46;47;48;49;50;51;52;55;56;57;58;59;60;61;62;69;71;72;73;75;76;100;101;102;103;104;105;106;107;108;109;112;113;114;115;116;117;118;119;120;121;122;123;124;125;128;129;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;151;152;153;154;155;156;157;160;161;162;163;164;165;166;167;168;169;170;171;172;173;174;175;176;177;178;179;180;181;182;183;184;185;186;187;190;191;192;193;194;195;196;197;198;199;200;201;202;203;204;205;206;207]';
Params.BrainRegions.NotGMorWM=[0 4 11 46 49 50 51 52 61 62];
Params.BrainRegions.WMregions=[35 40 41 44 45];
Params.BrainRegions.GMregions=AtlasVals;Params.BrainRegions.GMregions(ismember(Params.BrainRegions.GMregions,Params.BrainRegions.WMregions)...
    |ismember(Params.BrainRegions.GMregions,Params.BrainRegions.NotGMorWM))=[];
Params.BrainRegions.B0regions=[102 103 104 105 120 121 122 123 124 125 132 133 140 141 146 147 154 155 178 179 186 187 202 203 ...%OFC
31 32 47 48 116 117 170 171]; % amygdala/temporal lobe/hippocampus

end