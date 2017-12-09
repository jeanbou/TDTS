%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realize by default T-DTD initalizing when main while starting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- LIST OF THE CONSTANTS INDEPENDENT OF GUI --------------------- 
m_CONST_MAX_PRC = 100;                          % can be = 1, depends on the type of interpretation of maximal possible percentage ratio
M_CONST_CLASS_LABEL_UNDETERM = 777;             % label vector as class 777 in case when vector is not recognized
min_Size = 0.01;                                % minimal threshold step / to build optimal subthreshold of complexity interval in order to find optimal one for minimum of optimazing interations
zeroMeaningValue = 1e-012;                      % epsilon neighborhood
M_SP_COEF = 0.61803398875;                      % golden ratio / spliting ration reguired for optmila threshold search, can be 0.5 as well
M_CONS_MIN_THRESHOLD = 0.1;                     % ALFA for 1 - ALFA / to get sigma or subthreshold of PDF of Complexities of maximal possible database decomposition
m_RBFParams4ECTD.dstn_type = 'EUCL';            % Parameters required ZISC complexity estimator: 1) Type of distance
m_RBFParams4ECTD.MIF = 0;                       % 2) maximal influence field
m_RBFParams4ECTD.n_iteration = 1;               % 3) numbers of iteration applied in order to get average ratio of ZISC classification complexity
m_RBFParams4ECTD.m_pol_pow = 4;                 % 4) power of the polynome that approximates behaviour of the tail of the complexity ration
resolution_B = 2;                               % resolution required for PRISM methods based complexity estimator
M_CONST_OPTIMUM_STRFILENAME = 'm_optimum.mat';  % file which accumulates values controling self optimazing T-DTS
%%%%%%% GLOBAL VARIABLE THAT CONSISTS RESULT ANT ITS STATISTICS
M_CONST_STR_ALLRESULTSFILENAME  = 'm_results';
M_CONST_STR_STATRESULTSFILENAME = 'm_results_stat';
%%%%%%%% SEPARATES THREHOLD INTERVAL ON THE SUB INTERVAL TO BUILD A PDF OF COMPLEXITY
M_CONST_THRESHOLD_DELTA = 128-1;

%---------------- IN THIS PART WE INITIALIZED CONTROLS WITH SETTINGS BY DEFAULT
set(handles.figure1,'Name','T-DTS 2.50');
fprintf('\nT_DTS.m :: GUI started\n');
%%%% DU
set(findobj('Tag','m_du_menu'),'Value',1);
m_DU.method = 'CNN';
set(findobj('Tag','m_du_param_menu'),'Value',1);
m_DU.Dec_Param = 'Min_dist_to_prototype';
set(findobj('Tag','m_pu_menu'),'Value',1);
%%%% PU
m_PU.NoN = m_getDecMethod_params(m_DU.method);
m_PU.method = 'LVQ1';
%%%% SHOW PROCESSING
set(findobj('Tag','m_dsp_prcss_ch_bx'),'Value',0);             
m_Display_PARAMS.show_flag = 0;
m_Display_PARAMS.show_learningSlaves_batchSize = 2000;
m_detalized_PU_method = m_getProcMethod_params('none',0,0);
m_Display_PARAMS.show_learningSlaves_batchSize = m_detalized_PU_method.Epo; %%%% BY DEFAULT PARAMETERS CAN BE CHANGED
clear('m_detalized_PU_method');
%%%% RUN MODE SELECTION 
set(findobj('Tag','m_dsp_prcss_ch_bx'),'Value',1);
m_str_TDTS_mode = 'Uni';
%%%% NUMBER OF ITERATIONS
set(findobj('Tag','m_it_num_txt'),'String','1');
m_iterationNumber = 1;
%%%% COMPLEXITY THRESHOLD
set(findobj('Tag','m_cmplx_thr'),'String','0.500');
m_thresholdInterval = [0.5];
%%%% COMPLEXITY SELECTION
set(findobj('Tag','m_cmplx_lst'),'Value',1);
m_complexityESTList = {'Maximum_Standard_Deviation'};
%%%% INPUT FILENAME
set(findobj('Tag','db_filename'),'String','SQ2_2'); %%% by default
m_str_dataFileName = 'SQ2_2';
%%%% OUTPUT FILENAME
set(findobj('Tag','m_rslt_file'),'String','m_output'); %%% by default
m_str_outputFileName = 'm_output';
%%%% RANDOMIZATION FLAG
set(findobj('Tag','m_ran_data_ch_bx'),'Value',0);             
m_randomizeIncomingData_flag = 0;
%%%% LEARNING DATABASE EXTRACTION MODE SELECTION
set(findobj('Tag','m_lrn_db_extr_mode'),'Value',1);
m_LearningDBCreating.mode = 'Balanced';
%%%% SIZE OF LEARNING DATABASE in % SELECTED
set(findobj('Tag','m_lrn_db_rate'),'String','50.00');
m_LearningDBCreating.percentage = 50;
%%%% INIT NORMALIZATION / DATA CONVERTING MODE
set(findobj('Tag','m_data_conv'),'Value',1);
m_str_Type_ofNormalizing = 'No_Normalizing';
%%%% PCA/T
set(findobj('Tag','m_pca_ch_bx'),'Value',0);
m_PCA.flag = 0;
set(findobj('Tag','m_pca_prc'),'String','0.01');
m_PCA.per_elmnt_cmps = 1.0;
%%%% Initialize distance type selection
set(findobj('Tag','m_dist'),'Value',1);
m_str_Type_ofDistance = 'Euclidean';
save(M_STR_PARAMSFILENAME);
fprintf('T_DTS.m :: T-DTS initialized\n');