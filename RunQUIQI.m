function RunQUIQI(Subregions,Lambda,DataType,AnalType)
% Main function enabling the insertion of a motion degradation index into
% the analysis of MRI data. 
%
% INPUTS: 
%       - Subregions: cell array with the list of tissue types to be analysed 
%                       p1: grey matter
%                       p2: white matter
%       - Lambda: cell array with the list of the desired powers of the 
%                motion degradation index used for the computation of the 
%                ReML basis functions.
%                       0 : ordinary least-square analysis
%                       [0 1 2 3 4] : separate basis functions will be 
%                           computed from the MDI, one for each value.
%       - DataType: cell array with the list of maps to be analysed. 
%                   The options are:
%                       PDT1_R2s: R2* maps computed from the PDw and T1w raw image contrasts
%                       R2s_OLS: R2* maps computed from the PDw, T1w and MTw raw image contrasts
%                       R1: R1 maps computed from the T1w and PDw raw image contrasts
%                       MT: Mangetization Transfer maps computed from the PDw, T1w and MTw raw image contrasts
%       - AnalType: string with the name of the analysis to be performed. 
%                  The analysis runs on all the datasets in Subject_details.mat 
%                       Full - one-sample t-test of the age-dependence of the data on age
%                       Residuals - similar to 'Full' but here residual images 
%                           are saved for subsequent analysis
%                       Exclusion - conducts analysis after exclusion of the fraction
%                            of the cohort with the highest value of the MDI
%                       MotionBias - examines the dependence of the R2* data on the MDI (bias)
%                       Specificity - similar to 'Full' but the age assigned 
%                           to each dataset is scrambled randomly across the 
%                           datasets to allow for the monitoring of false positives
%
% OUTPUTS:
%       none
%
% EXAMPLES:
%
%       RunQUIQI({'p1'},{0},{'R1'},'Residuals')  
%       RunQUIQI({'p1'},{0},{'MT'},'Exclusion')  
%       RunQUIQI({'p1','p2'},{0},{'PDT1_R2s','R2s_OLS','R1','MT'},'Exclusion')  
%       RunQUIQI({'p1','p2'},{0,[0 1 2 3 4]},{'PDT1_R2s','R2s_OLS','R1','MT'},'Residuals') 
%       RunQUIQI({'p1','p2'},{0},{'PDT1_R2s','R2s_OLS','R1','MT'},'Exclusion') 
%       RunQUIQI({'p1','p2'},{0,[0 1 2 3 4]},{'PDT1_R2s','R2s_OLS','R1','MT'},'MotionBias') 
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland
% 
% REFERENCES
% 1. 'QUIQI – using a QUality Index for the analysis of Quantitative Imaging
% data'. Di Domenicantonio et al., Proceedings of the Annual Meeting of the
% International Society for Magnetic Resonance in Medicine, 2020.
%  2. 'Restoring statistical validity in group analyses of motion-corrupted MRI data',
% Lutti et al, Human Brain Mapping 2022.

% Get subject details
Params=GetParams;
eval(['load ' fullfile(Params.HomeDir,'Subject_Details.mat')]);

% Loop for all the DataType entries
for ctr=1:size(DataType,2)
    
    % Prepare analysis folders
    [QUIQI,~]=PrepAnalysis(Subject_Details,Lambda,DataType(ctr),Subregions,AnalType);
    
    % Run analysis
    RunAnalysis(QUIQI,AnalType);
    
    % Compare the MDI and the residuals analysis 
    MDIvsResAnalysis(QUIQI);
    
end


end