function [P,Rsq,yfit]=myPolyFit(X,Y,Powers,FitMethod)
% Polynomial fitting routine
%
% INPUTS:
%     - X: independent variable
%     - Y: dependent variable
%     - Powers: order of the polynomial fit
%     - FitMethod: allowed values: 1. 'NonNeg' - enforces positive
%     polynomial coefficients. 2. 'Free' - allows positive and negative
%     polynomial coefficients.

% OUTPUTS:
%     - P: polynomial coefficients
%     - Rsq: R-square of the fit
%     - yfit: fitted dependent variable
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Powers=linspace(0,Powers,Powers+1);
% P=zeros(size(Powers,2),size(Y,2));
Rsq=zeros(1,size(Y,2));
for ctr=1:size(Y,2)
    Ytemp=squeeze(Y(:,ctr));
    Xmat=[];
    for ctr2=1:size(Powers,2)
        if Powers(ctr2)==0
            Xmat=cat(2,Xmat,ones(size(X,1),1));
        else
            Xmat=cat(2,Xmat,X.^Powers(ctr2));
        end
    end
    if strcmp(FitMethod,'NonNeg')
        P{ctr} = lsqnonneg(Xmat,Ytemp);
    elseif strcmp(FitMethod,'Free')
        P{ctr}=pinv(Xmat'*Xmat)*Xmat'*Ytemp;
%         P(:,ctr)=polyfit(X,Ytemp,max(Powers));
%         yfit =polyval(P(:,ctr),X);%yfit=P(1)*R2sArray+P(2)
%         P(:,ctr)=P(end:-1:1,ctr)';%consistently with 'NonNeg', P is sorted in increasing powers of the model
    end
    yfit =Xmat* P{ctr}; %polyval(P(:,ctr),X);%yfit=P(1)*R2sArray+P(2)
    SSresid = sum((Ytemp - yfit).^2);
    SStotal = (length(Ytemp)-1) * var(Ytemp);
    Rsq(ctr) = 1 - SSresid/SStotal;%     Slope = P(1)
end

end