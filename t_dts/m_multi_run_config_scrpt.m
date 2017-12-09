%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script consists the set of parameters required for multi run mode of T-DTS        %
% This script calls when appropriate mode is selected                                    %
% PLESE MODIFY 3 MAIN PARAMETERS AND THEN PRESS SAVE, after select 'Multi' & press 'Run' %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% INIT COMPLEXITY ESTIMATOR LIST
% - m_complexityESTList = {'Maximum_Standard_Deviation','Fisher_Discriminant_Ratio','Purity_Measure','Normalized_mean_Distance','Divergence_Measure','Jeffries_Matusita_Distance','Bhattacharyya_Bound','Mahalanobis_Distance','Interclass_DM_CRT_Trace','Interclass_DM_CRT_Div_Trace','Interclass_DM_CRT_Div_Trace','Interclass_DM_CRT_Log_Det','Interclass_DM_CRT_DifTR','RBF_ZISC_based_Fusion_CRT','k_Nearest_Neighbors_Estimator','Collective_Entropy'};
m_complexityESTList = {'Fisher_Discriminant_Ratio'};
%%%% INIT THRESHOLDS INTERVAL
m_thresholdInterval = [0.0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
%%%% INIT NUMBER OF ITTERATTION
m_iterationNumber = 2; 