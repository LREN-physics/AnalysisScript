function MDIvsResAnalysis(QUIQI)
% Analysis of image noise vs motion degradation index. 
% For a given tissue type, compare the MDI values with the spatial variance 
% of the residuals obtained from the model under testing. Two measures of
% heteroscedasticity on the residuals are computed: a global (tissue type)
% and local (voxel level).
%
% INPUTS:
%     - QUIQI: structure containing all information used for analysis. Computed in PrepAnalysis.m.
%
% OUTPUTS:
%       none
%
% If paralell processing is prefered: uncomment lines 50-54 and comment line 56
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

MDIVals4Anal=cell(1,size(QUIQI,2));
Params=GetParams;
NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));

% Loop for the several analysis to be conducted. Such information is
% contained in QUIQI structure
for datactr=1:size(QUIQI,2)  
    
    eval(['load ' fullfile(QUIQI(datactr).CohortPath,'Subject_Details.mat')]);
    MDIVals=[];
    for ctr=1:size(Subject_Details,2)
        MDIVals=cat(1,MDIVals,[Subject_Details(ctr).QA.SDR2s.MTw Subject_Details(ctr).QA.SDR2s.PDw Subject_Details(ctr).QA.SDR2s.T1w]);        
    end

    MDIVals4Anal{datactr}=MDIVals(:,QUIQI(datactr).SDR2sIndx);
    
    CurrentPath=fullfile(QUIQI(datactr).CohortPath,QUIQI(datactr).AnalDir);
    P=spm_select('FPList',CurrentPath,'^ExplicitMask_.*.nii$');
    Vsave=spm_vol(P);Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ResMask.nii');
    ResMask=spm_read_vols(spm_vol(P));
    if contains(QUIQI(datactr).InputData,'PDw_R2s') && strcmp(QUIQI(datactr).TissueType,'p1')...
            || contains(QUIQI(datactr).InputData,'T1w_R2s') && strcmp(QUIQI(datactr).TissueType,'p1')...
            || contains(QUIQI(datactr).InputData,'MTw_R2s') && strcmp(QUIQI(datactr).TissueType,'p1')...
            || contains(QUIQI(datactr).InputData,'PDT1_R2s') && strcmp(QUIQI(datactr).TissueType,'p1')...
            || contains(QUIQI(datactr).InputData,'R2s_OLS') && strcmp(QUIQI(datactr).TissueType,'p1')
        ResMask(ismember(NMatlas,Params.BrainRegions.B0regions))=0;
    end
    spm_write_vol(Vsave,ResMask);
end

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool;
parfor datactr=1:size(QUIQI,2)

%for datactr=1:size(QUIQI,2)
    CurrentPath=fullfile(QUIQI(datactr).CohortPath,QUIQI(datactr).AnalDir);
    ResFiles=spm_select('FPList',CurrentPath,'^Res_.*.nii$');
    Vsave=spm_vol(ResFiles(1,:));
    
    ResMask=spm_read_vols(spm_vol(spm_select('FPList',CurrentPath,'^ResMask.nii$')));% for MDIvsRes analysis
    ResMaskIndx=find(ResMask~=0); 

    ResidVar=zeros(size(MDIVals4Anal{datactr},1),1);Residuals=zeros(size(ResMaskIndx,1),size(MDIVals4Anal{datactr},1));
    for Subjctr=1:size(MDIVals4Anal{datactr},1)% reads-in individual residual maps and estimates variance
        Subjctr
        tempRes=spm_read_vols(spm_vol(ResFiles(Subjctr,:)));
        Residuals(:,Subjctr)=tempRes(ResMaskIndx);% For voxel-wise heteroskedasticity analysis
        ResidVar(Subjctr)=var(Residuals(:,Subjctr),'omitnan');% For global heteroskedasticity analysis
    end
    SavePath=fullfile(spm_str_manip(Vsave.fname,'h'),'ResidualAnalysis');
    if ~exist(SavePath,'dir')
        mkdir(SavePath)
    end
    Vsave.fname=fullfile(SavePath,spm_str_manip(Vsave.fname,'t'));
    parsave(fullfile(SavePath,'ResidVar'),ResidVar)
    
    if ~isempty(strfind(QUIQI(datactr).AnalDir,'R2s'))%for relaxation rates in s-1
        ResidVar=ResidVar*1e6;
    elseif ~isempty(strfind(QUIQI(datactr).AnalDir,'R1'))
        ResidVar=ResidVar*1e-6;
    end
    ResFit=MDIvsResFit(SavePath,MDIVals4Anal{datactr},ResidVar);
    parsave(fullfile(SavePath,'MDIVals'),MDIVals4Anal{datactr})
    parsave(fullfile(SavePath,'ResFit'),ResFit)

    % KS & ARCH tests 
    % Only runs if Econometrics Toolbox exist
    
    if  license('test','Econometrics_Toolbox')~=0
        
        [~,B]=sort(ResFit,1);
        
        % Estimation of optimal lag from GLOBAL estimates of noise
        test=ResidVar(B);
        numLags = 4e1;
        logL = zeros(numLags,1); % Preallocation
        
        for k = 1:numLags
            Mdl = garch(0,k);                                    % Create ARCH(k) model
            [~,~,logL(k)] = estimate(Mdl,test,'Display','off');  % Obtain loglikelihood
        end
        
        aic = aicbic(logL,1:numLags);   % Get AIC
        [~,lags] = min(aic)             % Obtain suitable number of lags
        parsave(fullfile(SavePath,'lags'),lags)
        
        Residuals=Residuals';
        Hhet=zeros(size(Residuals,2),1);Phet=zeros(size(Residuals,2),1);        
        tic
        for ctr=1:size(Residuals,2)
            Signal=Residuals(:,ctr);
            if isempty(find(isnan(Signal)))&&isempty(find(isinf(Signal)))
                [Hhet(ctr),Phet(ctr)]=archtest(Signal(B)-mean(Signal),'Lags',lags);
            end
        end
        toc
        Hmap=zeros(size(spm_read_vols(spm_vol(ResFiles(1,:)))));
        Hmap(ResMaskIndx)=Hhet;
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'Hhet.nii');spm_write_vol(Vsave,Hmap);
        
        Pmap=zeros(size(spm_read_vols(spm_vol(ResFiles(1,:)))));
        Pmap(ResMaskIndx)=Phet;
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'Phet.nii');spm_write_vol(Vsave,Pmap);
        
    end
    
    % Delete residual individual residual maps to save disk space
    for ctr=1:size(ResFiles,1)
        delete(deblank(ResFiles(ctr,:)));
    end
        
end

end

function ResFit=MDIvsResFit(SavePath,MDI,Res)
Params=GetParams();

if ~exist(SavePath,'dir')
    mkdir(SavePath)
end

[P,Rsq,ResFit]=myPolyFit(MDI,Res,Params.MDIvsResOrder,'Free');
FittingEstimates.P=P;FittingEstimates.Rsq=Rsq;
save(fullfile(SavePath,'FitEstimates.mat'),'FittingEstimates', '-v7.3')

end
function parsave(Path,Var)
save(Path,'Var')
end


