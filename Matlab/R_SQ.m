function [R_SQ] = R_SQ(Predicted,Actual)
%R_SQ Computes the R_SQ value for a given set of predicted and actual data
%   INPUTS:
%         Predicted = Predicted Data
%         Actual = Actual Data
%   OUTPUTS:
%         R_SQ = R^2 correlation coefficient
SSR = sum((Predicted - Actual).^2);
TSS = sum(((Actual - mean(Actual)).^2));
R_SQ = 1 - SSR/TSS;
end