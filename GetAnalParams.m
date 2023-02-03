function AnalParams=GetAnalParams(AnalType)
% Retrieves parameters specific to each analysis type.
%
% INPUTS: 
%     - AnalType: string label of the analysis to be conducted.
% OUTPUTS:
%     - AnalParams: structure of analysis parameters. 
%
% Default analysis settings are set for a standard 'Full' analysis of the whole cohort.
% For other analysis types, the relevant parameters are set specifically.
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

% Set default analysis settings for a standard 'Full' analysis of the whole cohort.
% For other analysis types, the relevant parameters are set specifically.
AnalParams.SaveResiduals=0;         % 1 to save residual maps; 0 otherwise
AnalParams.Subcohort=0;             % 1 if one wants to use a subcohort of the available data; 0 otherwise
AnalParams.Exclusion=0;             % 1 if one wants to exclude some percentage of the data; 0 otherwise
AnalParams.AgeBinning.Nsamples=Inf; % Inf if we want to use all the data available; finite number if one wishes to use just a few samples of the whole cohort
AnalParams.NRepeats=1;              % number of analysis repetitions
AnalParams.AgeBinning.Agemin=10;    % min age for the binning of the data
AnalParams.AgeBinning.Agemax=100;   % max age for the binning of the data
AnalParams.AgeBinning.BinNb=19;     % number of bins for the binning of the data
AnalParams.ShuffleAge=0;            % 1 if one wishes to shuffle the age across the data; 0 otherwise
AnalParams.ShuffleData=0;           % 1 if one wishes to shuffle the data across subgroups for Group Comparison; 0 otherwise 
AnalParams.NGroup1=0;               % number of datasets in the first of the two groups created with the selected data for Group Comparison

% % To create public sample:
% AnalParams.SaveResiduals=0;AnalParams.Subcohort=1;AnalParams.Exclusion=0;
% AnalParams.AgeBinning.Nsamples=10;AnalParams.NRepeats=1;
% AnalParams.AgeBinning.Agemin=10;AnalParams.AgeBinning.Agemax=100;AnalParams.AgeBinning.BinNb=19;
% AnalParams.ShuffleAge=0;AnalParams.ShuffleData=0;AnalParams.NGroup1=0;

% Set analysis parameters depending on the type of analysis:
if strcmp(AnalType,'Exclusion')
    AnalParams.SaveResiduals=1;                         % save residual maps
    AnalParams.Exclusion=[0 3 7 13 20 30];              % percentage of the data to be excluded
    AnalParams.NRepeats=size(AnalParams.Exclusion,2);   % repeats the analysis for each percentage value defined above

elseif strcmp(AnalType,'Residuals')
    AnalParams.SaveResiduals=1;                         % save residual maps
    
elseif strcmp(AnalType,'Reproducibility')
    AnalParams.Subcohort=1;                             % use a subcohort of the data  
    AnalParams.AgeBinning.Nsamples=200;                 % number of datasets to be selected for the cohort
    AnalParams.NRepeats=100;                            % number of repetitions of the analysis
    
elseif strcmp(AnalType,'Specificity')
    AnalParams.AgeBinning.Nsamples=10;                  % number of datasets to be selected for the cohort
    AnalParams.NRepeats=1e3;                            % number of repetitions of the analysis
    AnalParams.Subcohort=1;                             % use a subcohort of the data  
    AnalParams.ShuffleAge=1;                            % shuffle the age between the datasets
    
elseif strcmp(AnalType,'MotionBias')
    AnalParams.SaveResiduals=1;
%     Commented out for the released data, which is already a subset of the
%     cohort used for publication
%     AnalParams.Subcohort=1;AnalParams.AgeBinning.Nsamples=200;
%     AnalParams.AgeBinning.Agemin=56;AnalParams.AgeBinning.Agemax=58;AnalParams.AgeBinning.BinNb=2;

elseif strcmp(AnalType,'GroupComparison')
    AnalParams.Subcohort=1;                             % use a subcohort of the data  
    AnalParams.ShuffleData=1;                           % shuffles the data between the subgroups of subjects
    AnalParams.NGroup1=10;                              % number of datasets in the first of the two groups created with the selected data 
    AnalParams.AgeBinning.Nsamples=200;                 % number of datasets to be selected for the cohort 
    AnalParams.NRepeats=1e3;                            % number of repetitions of the analysis
    AnalParams.AgeBinning.Agemin=56;                    % min age for the subcohort
    AnalParams.AgeBinning.Agemax=58;                    % max age for the subcohort
    AnalParams.AgeBinning.BinNb=2;                      % number of samples per bin
end
end