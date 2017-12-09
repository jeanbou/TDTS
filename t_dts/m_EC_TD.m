function [di,cm] = m_EC_TD(X,T,numPrototypes,m_Complexity_PARAMS,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);
%----------------------------------------- dMLDiv.m :: function DMLDiv(X,T) ----------------------------------------------%
% This function estimate complexity of the problem X, labed by classes T, than obtained complexity ratio, send to decision%
% function which returns flag of decision; 1 - means to devide, 0 - means stop it                                         %
% m_complexity_PARAMS.Consists - Threshold value of complexity                                                            %                                                                     %
% zeroMeaningValue - ration that's closed to 0, but not equal, requires for comparizing float output to zero              % 
% m_RBFParams4ECTD - structure of two parameters, which controls mZISC complexity estimating process                      %
% resolution_B _ parameter which responce for the resolution of Purity and Collective_Entropy procedures, by PRISM Singh  %
% The function outputs :                                                                                                  %
%   * divide decision parameter 'di'                                                                                      %
%       ** if 'di' equals 1 it means to devide                                                                            %
%       ** if 'di' equals 0 stop deviding process                                                                         %
%   * complexity value parameter 'cm'                                                                                     %
%       ** real interval [0, 1]                                                                                           %
% It calls directly the complexity estimation functions                                                                   %
%    which are adoped internally to multi-class cases                                                                     %
%-------------------------------------------------------------------------------------------------------------------------%

if ( max(T) == min(T) )
    if m_Display_PARAMS.show_flag
	    fprintf('m_EC_D.m :: Note(1) :: Pure cluster detected\n');
    end;
    di = 0; %%% Stop decomposition
    cm = 1; %%% Predefine problem as simple 1(non)-complex
    return;
end;

if ( size(T,2) <= numPrototypes )
    if m_Display_PARAMS.show_flag
	    fprintf('m_EC_D.m :: Note(2) :: Minimal cluster size detected\n');
    end;
    di = 0; %%% Stop decomposition
    cm = 1; %%% Predefine problem as simple 1(non)-complex
    return;
end;

switch m_Complexity_PARAMS.Method_name,
	case 'Maximum_Standard_Deviation',   cm = m_countMaxstd(X,T);
    case 'Fisher_Discriminant_Ratio',    cm = m_FisherDR(X,T);  
    case 'Purity_Measure',               cm = m_PurityM(X,T,resolution_B);
    case 'Normalized_mean_Distance',     cm = m_NormalizedMD(X,T); 
    case 'KLD',                          cm = m_KLD(X,T); 
    case 'JMDBC',                        cm = m_JMDBC(X,T);
    case 'Bhattacharyya_Coefficient',    cm = m_BhattacharyyaCoefficient(X,T);
    case 'Mahalanobis_Distance',         cm = m_MahalanobisD(X,T,zeroMeaningValue);
    case 'Interclass_DM_CRT_Trace',      cm = m_IDM_CR_TrINV(X,T,zeroMeaningValue);
    case 'Interclass_DM_CRT_Div_Trace',  cm = m_IDM_CR_DivTR(X,T,zeroMeaningValue);
    case 'Interclass_DM_CRT_Log_Det',    cm = m_IDM_CR_LogDetINV(X,T,zeroMeaningValue);
    case 'Interclass_DM_CRT_DifTR',      cm = m_IDM_CR_DifTR(X,T);
    case 'RBF_ZISC_based_Fusion_CRT',    cm = m_RBF_ZISC_fusion(X,T,m_RBFParams4ECTD.dstn_type,m_RBFParams4ECTD.MIF,m_RBFParams4ECTD.n_iteration,m_RBFParams4ECTD.m_pol_pow,m_Display_PARAMS);
    case 'k_Nearest_Neighbors_Estimator',cm = m_kNN_Estimator(X,T,resolution_B);
    case 'Collective_Entropy',           cm = m_Collective_Entropy(X,T,resolution_B);
    case 'JSD',                          cm = m_JSD(X,T);
    case 'Hellinger_Distance',           cm = m_Hellinger_Distance(X,T);
    case 'Bayes_Error',                  cm = m_bayes_error(X,T);
        
    otherwise
        cm =NaN
        fprintf('m_EC_TD.m :: Error(1) :: No method description for selected complexity estimator\n');
	    return;
end;

%%%%%%%%%%%%%%% DATA FILTER :: Check complexity estimator %%%%%%%%%%%%%%
if (cm > 1+zeroMeaningValue) | (cm < -zeroMeaningValue)
    fprintf('m_EC_TD.m :: Error(2) :: Complexity measure value %2.3f is out of range [0;1]\n',cm);
    return;
end;

%%%%%%%%%%%%% Taking decision :: Simple procedure, can be done more sofisticated %%%%%%%%%%%%
if (m_Complexity_PARAMS.Threshold > cm) %%%% fixed threshold ration is higher (less complex) than it is required
    di = 1;
    if m_Display_PARAMS.show_flag
        fprintf('m_EC_TD.m :: Complexity: %2.3f < Threshold: %2.3f\n', cm, m_Complexity_PARAMS.Threshold);
    end;
else
	di = 0; %%% Stop Dividing process
end;