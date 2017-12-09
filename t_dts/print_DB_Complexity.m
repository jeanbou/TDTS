function print_DB_Complexity(m_str_dataFileName,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS, m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
%------------------------------------------------------------------------------------------------------------%
% This function prints out all avaliable complexities for solid /not partioned/ database                     %
% which set by default. Linked to GUI                                                                        %
% minSize - the size of X, which cannot be decomposed                                                        %
% zeroMeaningValue - ration that's closed to 0, but not equal, requires for comparizing float output to zero % 
% m_RBFParams4ECTD - structure of two parameters, which controls mZISC complexity estimating process         %
%------------------------------------------------------------------------------------------------------------%

[inPL, inCL, inPG, inCG] = m_uploadDBs(m_str_dataFileName, m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
X = [inPL, inPG];
T = [inCL, inCG];
clear('inPL', 'inPG','inCL','inCG');

cm = m_countMaxstd(X,T);
fprintf('m_print_DB_Coplexity :: Maximum of the Standard Deviation :: %1.6f\n',cm);

cm = m_FisherDR(X,T);
fprintf('m_print_DB_Coplexity :: Fisher Discriminant Ratio :: %1.6f\n',cm);

cm = m_PurityM(X,T,resolution_B);
fprintf('m_print_DB_Coplexity :: Purity Measure :: %1.6f\n',cm);

cm = m_NormalizedMD(X,T);
fprintf('m_print_DB_Coplexity :: Normalized Mean Distance :: %1.6f\n',cm);

cm = m_KLD(X,T);
fprintf('m_print_DB_Coplexity :: Kullback Divergence Measure :: %1.6f\n',cm);

cm = m_JMDBC(X,T);
fprintf('m_print_DB_Coplexity :: Jeffries-Matusita Distance Based Criterion :: %1.6f\n',cm);

cm = m_BhattacharyyaCoefficient(X,T);
fprintf('m_print_DB_Coplexity :: Bhattacharyya Coefficient :: %1.6f\n',cm);

cm = m_MahalanobisD(X,T,zeroMeaningValue);
fprintf('m_print_DB_Coplexity :: Mahalanobis_Distance :: %1.6f\n',cm);

cm = m_IDM_CR_TrINV(X,T,zeroMeaningValue);
fprintf('m_print_DB_Coplexity :: Interclass_DM_CRT_Trace :: %1.6f\n',cm);

cm = m_IDM_CR_DivTR(X,T,zeroMeaningValue);
fprintf('m_print_DB_Coplexity :: Interclass_DM_CRT_Div_Trace :: %1.6f\n',cm);

cm = m_IDM_CR_LogDetINV(X,T,zeroMeaningValue);
fprintf('m_print_DB_Coplexity :: Interclass_DM_CRT_Log_Det :: %1.6f\n',cm);

cm = m_IDM_CR_DifTR(X,T);
fprintf('m_print_DB_Coplexity :: Interclass_DM_CRT_DifTR :: %1.6f\n',cm);

cm = m_RBF_ZISC_fusion(X,T,m_RBFParams4ECTD.dstn_type,m_RBFParams4ECTD.MIF,m_RBFParams4ECTD.n_iteration,m_RBFParams4ECTD.m_pol_pow,m_Display_PARAMS);
fprintf('m_print_DB_Coplexity :: RBF_ZISC_based_Fusion_CRT :: %1.6f\n',cm);

cm = m_kNN_Estimator(X,T,resolution_B);
fprintf('m_print_DB_Coplexity :: k_Nearest_Neighbors_Estimator :: %1.6f\n',cm);

cm = m_Collective_Entropy(X,T,resolution_B);
fprintf('m_print_DB_Coplexity :: Collective_Entropy :: %1.6f\n',cm);

cm = m_JSD(X,T);
fprintf('m_print_DB_Coplexity :: Jensen-Shannon Divergence :: %1.6f\n',cm);

cm = m_Hellinger_Distance(X,T);
fprintf('m_print_DB_Coplexity :: Hellinger Distance :: %1.6f\n',cm);

cm = m_bayes_error(X,T);
fprintf('m_print_DB_Coplexity :: Bayes Error :: %1.6f\n',cm);