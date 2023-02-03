function [QUIQI,Folders]=PrepAnalysis(Subject_Details_All,lambda,DataTypes,Subregions,AnalType)
% Compiles into QUIQI and Folders structures all information required for
% subsequent analysis. 
% Also, edits cohort details if necessary (e.g. for Exclusion or Specificity analyses)
% and computes an explicit mask for image analysis if it doesn't exist (PrepareMasksAndCohort function)
%
% INPUTS: 
%     - Subject_Details_All: subject demographic information structure.
%     - DataTypes: type of input data for the analysis. Initialized in RunQUIQI.m.
%     - Subregions: region of interest of the analysis. Initialized in RunQUIQI.m.
%     - AnalType: analysis type. Initialized in RunQUIQI.m.
%
% OUTPUTS: 
%       QUIQI.TissueType: tissue type of interest in the analysis
%       QUIQI.SDR2sIndx: index of the MDI of interest in the demographic information structure, which
% contains three entries. Matches the type of data to be analysed i.e. 1/2/3 for MT/PD/T1-weighted data
%       QUIQI.MotionReg: MDI values extracted from the demographic information structure. From the 
% three entries, the extracted values are taken from the value of QUIQI.SDR2sIndx. Used for the MotionBias analyses.
%       QUIQI.ReML: set of basis functions for ReML
%       QUIQI.ROI: region of interest
%       QUIQI.AnalDir: analysis output folder
%       QUIQI.InputData: input data string 
%
% If paralell processing is prefered: uncomment lines 113-117 and comment line 119
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params=GetParams;                        % data paths and general parameters  
AnalParams=GetAnalParams(AnalType);      % paramaters specific to the anlysis type
RegionStr=RegionLabeltoName(Subregions); % name of the regions of interest for the analysis
NbRepeats=AnalParams.NRepeats;           % number of repeats for the analysis

RootPath=fullfile(Params.HomeDir,[Params.AnalysisFolder '_' AnalType]);
if ~exist(RootPath,'dir')
    mkdir(RootPath)
end
save(fullfile(RootPath,'AnalParams'),'AnalParams', '-v7.3')

%% Cohort definition

CreateNewMask=zeros(NbRepeats,1);CohortPath=cell(1,NbRepeats);

for subsetctr=1:NbRepeats % loops for repetitions
    CohortPath{subsetctr}=fullfile(RootPath,num2str(subsetctr));
    if ~exist(CohortPath{subsetctr},'dir')
        mkdir(CohortPath{subsetctr})
    end
    
    if isempty(spm_select('FPList',CohortPath{subsetctr},'^Subject_Details.*.mat$')) && ~strcmp(AnalType,'Exclusion') ...
            || strcmp(AnalType,'Exclusion') % With Exclusion analyses, the subject detail list is updated systematically as it varies with the type of data to be analysed
       
        if ~isinf(AnalParams.AgeBinning.Nsamples)
            Subject_Details=Select_cohort_subset(Subject_Details_All,AnalParams);
% %             Hack to enable selection of the same subset of data in	
% %             exclusion analyses ('Subject_Details_subset_for_Exclusion.mat')	
%             eval(['load ' fullfile(Params.DataDir,'Subject_Details_subset_for_Exclusion.mat')])	
%             Subject_Details = Subject_Details_bkup;
        else
            Subject_Details=Subject_Details_All;
        end
        
        if strcmp(AnalType,'Exclusion') 
            if size(DataTypes,2)>1
                error("Exclusion analyses should be conducted with each data type consecutively")
            elseif AnalParams.Exclusion(subsetctr)~=0
                MDIVals=zeros(size(Subject_Details,2),1);
                for subjctr=1:size(Subject_Details,2)
                    if strcmp(DataTypes{:},'PDw_R2s')
                        MDIVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.PDw;
                    elseif strcmp(DataTypes{:},'T1w_R2s')
                        MDIVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.T1w;
                    elseif strcmp(DataTypes{:},'MTw_R2s')
                        MDIVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.MTw;
                    elseif strcmp(DataTypes{:},'PDT1_R2s') ||  strcmp(DataTypes{:},'R1')
                        MDIVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.PDw + Subject_Details(subjctr).QA.SDR2s.T1w;
                    elseif strcmp(DataTypes{:},'R2s_OLS') ||  strcmp(DataTypes{:},'MT')
                        MDIVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.MTw + Subject_Details(subjctr).QA.SDR2s.PDw + Subject_Details(subjctr).QA.SDR2s.T1w;
                    end
                end
                [~,B]=sort(MDIVals);
                Indx=B(end-round(size(MDIVals,1)*AnalParams.Exclusion(subsetctr)/100)+1:end);
                Subject_Details(Indx)=[];
            end
        end
        
        if AnalParams.ShuffleAge
% %             Random permutation of the subject's age	
%             Age=[Subject_Details(:).Age];	
%             ShuffledAge=Age(randperm(size(Subject_Details,2)));	
%             for ctr=1:size(Subject_Details,2)	
%                 Subject_Details(ctr).Age=ShuffledAge(ctr);	
%             end

%           Subject's age is randomly assigned from a uniform distribution
%           ranging from the min and max age in the whole cohort.
            RandomAge=min([Subject_Details_All.Age])+(max([Subject_Details_All.Age])-min([Subject_Details_All.Age]))*rand(1,size(Subject_Details,2));
            for ctr=1:size(Subject_Details,2)
                Subject_Details(ctr).Age=RandomAge(ctr);
            end
            
        elseif AnalParams.ShuffleData
            Subject_Details=Subject_Details(randperm(size(Subject_Details,2)));
        end
        CreateNewMask(subsetctr)=1;
        save(fullfile(CohortPath{subsetctr},'Subject_Details'),'Subject_Details', '-v7.3')
    end
end

%% Mask creation

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool(12);
parfor subsetctr=1:NbRepeats

%for subsetctr=1:NbRepeats % loops for repetitions
    PrepareMasksAndCohort(CohortPath{subsetctr},CreateNewMask(subsetctr),Subregions)
end

Folders.CohortPaths=CohortPath;
Folders.AnalFolders={};Folders.DataFolders={};

%% Build QUIQI structure

ctr1d=0;
for subsetctr=1:NbRepeats % loops for repetitions
    eval(['load ' spm_select('FPList',CohortPath{subsetctr},'^Subject_Details.*.mat$')]);
    SDR2sVals=zeros(size(Subject_Details,2),3);
    for subjctr=1:size(Subject_Details,2)
        SDR2sVals(subjctr,:)=[Subject_Details(subjctr).QA.SDR2s.MTw Subject_Details(subjctr).QA.SDR2s.PDw Subject_Details(subjctr).QA.SDR2s.T1w];
    end
    for ctr1=1:size(DataTypes,2)
        for ctr=1:size(Subregions,2)
            for ctr2=1:size(lambda,2)
                ctr1d=ctr1d+1;
                QUIQI(ctr1d).CohortPath=CohortPath{subsetctr};
                if ischar(Subregions{ctr})
                    QUIQI(ctr1d).TissueType=Subregions(ctr);
                elseif isnumeric(Subregions{ctr})
                    if ismember(Subregions{ctr},Params.BrainRegions.GMregions)==1
                        QUIQI(ctr1d).TissueType={'p1'};
                    elseif ismember(Subregions{ctr},Params.BrainRegions.WMregions)==1
                        QUIQI(ctr1d).TissueType={'p2'};
                    end
                end
                if strcmp(DataTypes{ctr1},'MTw_R2s')
                    QUIQI(ctr1d).SDR2sIndx=1;
                elseif strcmp(DataTypes{ctr1},'PDw_R2s')
                    QUIQI(ctr1d).SDR2sIndx=2;
                elseif strcmp(DataTypes{ctr1},'T1w_R2s')
                    QUIQI(ctr1d).SDR2sIndx=3;
                elseif strcmp(DataTypes{ctr1},'PDT1_R2s')
                    QUIQI(ctr1d).SDR2sIndx=[2 3];
                elseif strcmp(DataTypes{ctr1},'R2s_OLS')
                    QUIQI(ctr1d).SDR2sIndx=[1 2 3];
                elseif strcmp(DataTypes{ctr1},'R1')
                    QUIQI(ctr1d).SDR2sIndx=[2 3];
                elseif strcmp(DataTypes{ctr1},'MT')
                    QUIQI(ctr1d).SDR2sIndx=[1 2 3];
                end
                QUIQI(ctr1d).MotionReg={SDR2sVals(:,QUIQI(ctr1d).SDR2sIndx)};
                if isnumeric(lambda{ctr2})
                    TempReML=[];
                    for ctr3=1:size(lambda{ctr2},2)
                        if lambda{ctr2}(ctr3)==0% to avoid duplicates of the identity matrix
                            TempReML=cat(2,TempReML,ones(size(SDR2sVals,1),1));
                        else
                            TempReML=cat(2,TempReML,SDR2sVals(:,QUIQI(ctr1d).SDR2sIndx).^lambda{ctr2}(ctr3));
                        end
                    end
                    QUIQI(ctr1d).ReML=TempReML;
                    LambdaStr=num2str(lambda{ctr2});LambdaStr=LambdaStr(find(isspace(LambdaStr)==0))';FolderName={};
                    for strctr=1:size(LambdaStr,1)
                        FolderName{strctr}=LambdaStr(strctr);
                    end
                    FolderName=['_lambda_' strjoin(FolderName,'_')];
                end
                QUIQI(ctr1d).ROI=RegionStr{ctr};
                QUIQI(ctr1d).AnalDir=fullfile(DataTypes{ctr1}, [char(RegionStr{ctr}) FolderName]);
                QUIQI(ctr1d).InputData=['fin_dart_' char(QUIQI(ctr1d).TissueType) '_ws.*' DataTypes{ctr1} '.nii$'];
                if ctr1==1&&subsetctr==1
                    Folders.AnalFolders(size(Folders.AnalFolders,2)+1)={[char(RegionStr{ctr}) FolderName]};
                end
            end
        end
    end
end
Folders.DataFolders=DataTypes;

end


