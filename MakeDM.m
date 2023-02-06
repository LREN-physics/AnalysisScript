function DM=MakeDM(Subject_Details,MDI,AnalType)
% Computes a design  matrix DM used in subsequent analyses (RunAnalysis.m).
%
% INPUTS:
%     - Subject_Details: structure of demographic information for the analysis cohort e.g. age,
% gender, brain volume,... 
%     - MDI: Motion Degradation Index values for each dataset of the analysis cohort. Inserted in the design matrix for
% 'MotionBias' analyses, for each power value of Params.MotionRegPowers 
%     - AnalType: string of the analysis type 

% OUTPUTS:
%     - DM.mat: design matrix data
%     - DM.size: size. Used for the definition of the F-contrast in the
%     ageing analyses 
%     - DM.text: text description
%     - DM.ReMLFcontrast: F-contrast used by ReML for identification of significant voxels
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland


Params=GetParams;
DM.confounds=[Subject_Details.confound];
if strcmp(AnalType,'MotionBias')
    X=MDI.^Params.MotionRegPowers(1);
    for ctr=2:size(Params.MotionRegPowers,2)
        X=cat(2,X,MDI.^Params.MotionRegPowers(ctr));    
    end
    DM.size=size(X,2);
    for ctr=1:size(X,2)
        DM.text{ctr}=num2str(ctr);
    end
else
    if Params.AgeModelOrder==3
        X=cat(2,[Subject_Details(:).Age]',[Subject_Details(:).Age].^2',[Subject_Details(:).Age].^3');
        DM.text={'Age','Agesq','AgeCub'};
    elseif Params.AgeModelOrder==2
        X=cat(2,[Subject_Details(:).Age]',[Subject_Details(:).Age].^2');%         DM.mat=cat(2,[Subject_Details(:).Age]',([Subject_Details(:).Age]-mean([Subject_Details.Age])).^2');
        DM.text={'Age','Agesq'};
    end
    DM.size=size(DM.text,2);
end
% Adds confounds to the design matrix
X=cat(2,X,[DM.confounds.Gender]',[DM.confounds.BrainVol]');
DM.text(DM.size+1)={'Gender'};DM.text(DM.size+2)={'BrainVol'};

% Orthogonalize design matrix
DM.mat=spm_orth(X-repmat(mean(X,1),[size(X,1) 1]));

DM.desc=AnalType;
end
