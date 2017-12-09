function m_TDTS(m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, M_CONST_THRESHOLD_DELTA,m_str_optThresholdFN,M_SP_COEF,M_CONST_ALFA,dname,m_Display_PARAMS,m_randomizeIncomingData_flag,M_CONST_SUBTHRESHOLD_STEP,zeroMeaningValue,m_RBFParams4ECTD,m_str_Type_ofDistance,m_str_Type_ofNormalizing,m_PCA,m_DU,m_PU,m_tmp_str_TDTS_mode,m_str_outputFileName,m_iterationNumber,m_thresholdInterval,m_complexityESTList,resolution_B,m_LearningDBCreating);
% ------------------------------------------------------------------------------------------------------------
% This function is main body-frame of the program, manages the run 3 main modes of the running of TDTS
% ------------------------------------------------------------------------------------------------------------
% This function consists sub-function m_itT_DTS(
% This function consists sub-function m_recDBSplit(
% This function consists sub-function m_recDBSplitFollowing_tree(

clear('m_results','m_results_stat');
m_results = [];
m_results_stat = [];

switch m_tmp_str_TDTS_mode
    
  case {'Uni','Multi'} % the only one difference between modes was that uploaded script of multi run configuration, with new range of parameters
  
    % If it's required complex multi run, than modify parameters and run it as 'UNI' mode, but with the list of new parameters that indicated in script
    clear('m_str_optThresholdFN');
    if strcmp(m_tmp_str_TDTS_mode,'Multi')
        m_multi_run_config_scrpt;
        if m_Display_PARAMS.show_flag
            fprintf('m_TDTS.m :: Parameters for Multi mode run has been picked up from script\n');
        end;
    end;
    for ind_selected_Complexity = 1 : size(m_complexityESTList,2)
        m_Complexity_PARAMS.Method_name = m_complexityESTList{ind_selected_Complexity};
        m_thrsh_index = 0; %%% Temporal index of Threshold from interval
        for flt_Threshold = m_thresholdInterval,  
            m_thrsh_index = m_thrsh_index + 1;
            m_Complexity_PARAMS.Threshold = flt_Threshold;
            %% if exist file where accumulated some results, so download this and then continue
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                clear('m_results', 'm_results_stat');
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;
            %%% launch ittertive T-DTS
            [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, COUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_thrsh_index, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);
            %%% saving results to last mat-file of the range of output result files
            save(strcat(m_str_outputFileName,'.mat'),'m_results', 'm_results_stat','POUT','COUT','WOUT','TREE','MAXTREELEVEL');
        end;
    end;
    
  case 'Self_Tuning',  
        
        ind_selected_Complexity = size(m_complexityESTList,2);   % we take last and in general case and the first complexity type from the lists // supose this list consist only one complexity estimator type
        m_thrsh_index = 1;
        m_Complexity_PARAMS.Method_name = m_complexityESTList{m_thrsh_index};
        if  (exist(m_str_optThresholdFN) == 2 ) %%% there is a file that consist parameters of self optimization so it means that we pick-up last result from it and continue to sear for optimum
            load(m_str_optThresholdFN);
            %% if exist file where accumulated some results, so download this and then continue utilizing info of self-optimizing and lasts results
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;            
        else 
            %% This is the first run of self optimization procedure // figure out the state of the task
            fprintf('m_TDTS.m :: Self-tuning :: Building optimal threshold sub interval\n');
            m_Complexity_PARAMS.Threshold = 1; %%% do first a maximal possible decomposition
            clear('m_genCML_List','m_cmpl_Index','m_numNonComplxCluster','m_treeLevel');
            m_genCML_List = [];
            %%% Obtaine complexity estimators list in order indicate complexity of the task and begin a wise search for optimal threshold
            clear('STR','str','WOUT','POUT','TOUT','m_CML_List','maxlevel');
            for i = 1 : m_iterationNumber
                STR = [];       % complex struction that accumulates information about decomposition
                m_CML_List = [];% list of the complexity thresholds appeared during decomposition 
                str='STR';      % additional helpful string that accumulates expression of decomposition
                WOUT = [];      % matrix of coordinates of the prototypes after decomposition
                POUT = [];      % structure of the decompossed data vectors
                TOUT = [];      % structure of the decompossed class labels 
                
                [inPL, inCL, inPG, inCG] = m_uploadDBs(dname, m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
                clear('inPG','inCG');
                [WN, POUT, TOUT, STR, maxlevel, m_CML_List] = m_recDBSplit(WOUT,POUT,TOUT,0,inPL,m_getDecMethod_params(m_DU.method),inCL,1,STR,str,m_str_Type_ofDistance,m_Complexity_PARAMS,m_CML_List,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);                     
                m_CML_List = sort(m_CML_List);
                m_cmpl_Index = 1;
                while ( m_CML_List(m_cmpl_Index) < 1 )
                    m_genCML_List = [m_genCML_List, m_CML_List(m_cmpl_Index)];
                    m_cmpl_Index = m_cmpl_Index + 1;
                end;
                m_numNonComplxCluster(i) = size(m_CML_List,2) - m_cmpl_Index + 1;
                m_treeLevel(i) = maxlevel;
                clear('STR','str','WOUT','POUT','TOUT','m_CML_List','maxlevel');
            end;            
            m_avgNumSimCluster = mean(m_numNonComplxCluster);
            m_stdNumSimCluster =  std(m_numNonComplxCluster);
            m_avgTreeLevel = mean(m_treeLevel);
            m_stdTreeLevel =  std(m_treeLevel);
            clear('m_numNonComplxCluster','m_treeLevel');
            m_genCML_List = sort(m_genCML_List);
            min_cmp_Threshold = min(m_genCML_List);
            max_cmp_Threshold = max(m_genCML_List);
            length_ComplexityInterval = max_cmp_Threshold - min_cmp_Threshold;
            if ( length_ComplexityInterval >= M_CONST_SUBTHRESHOLD_STEP )
                m_histXInterval = ( min_cmp_Threshold : (max_cmp_Threshold-min_cmp_Threshold)/M_CONST_THRESHOLD_DELTA : max_cmp_Threshold );
            else
                error('m_TDTS.m :: Error(1) :: length_ComplexityInterval is smaller than M_CONST_SUBTHRESHOLD_STEP');
            end;
            for i = 1 : ( size(m_histXInterval,2) - 1 )
                m_ind_histXInterval = find ( m_histXInterval(i) <= m_genCML_List & m_genCML_List < m_histXInterval(i+1) ); 
                m_pdfCML(i) = size(m_ind_histXInterval,2)/size(m_genCML_List,2);
            end;
            m_pdfCML(size(m_histXInterval,2)) = 1 - sum(m_pdfCML);
            clear('m_ind_histXInterval','i');
            %%%% Find sub-interwal of m_genCML_List where SUM(PDF) = 1 - ALFA
            m_SUMPDF = max(m_pdfCML);
            m_rightIndex = max(find(m_SUMPDF == m_pdfCML));
            m_leftIndex = m_rightIndex;
            m_leftIndexValue = m_SUMPDF;
            m_rightIndexValue = m_SUMPDF;
            m_reachedLeftBorderflag = 0;
            m_reachedRightBorderflag = 0;
            while ( m_SUMPDF < ( 1 - M_CONST_ALFA ) )
                %%%% Stop to go on left or right
                if (m_leftIndex == 1 )
                    m_reachedLeftBorderflag = 1;
                end;
                if (m_rightIndex == size(m_pdfCML,2) )
                    m_reachedRightBorderflag = 1;
                end;
                %%%% Look loock-up if you has not reached a border what is on the left or on the right
                if ~m_reachedLeftBorderflag
                    m_leftIndex = m_leftIndex - 1;
                    m_leftIndexValue = m_pdfCML(m_leftIndex);
                end;
                if ~m_reachedRightBorderflag
                    m_rightIndex = m_rightIndex + 1;
                    m_rightIndexValue = m_pdfCML(m_rightIndex);
                end;
                %%% Now move and increase sum
                if ( (m_leftIndexValue >= m_rightIndexValue ) & ~m_reachedLeftBorderflag )
                    m_SUMPDF = m_SUMPDF + m_leftIndexValue;
                    %%% we do not go on the right
                    if  ~m_reachedRightBorderflag
                        m_rightIndex = m_rightIndex - 1;
                    end;
                else
                    m_SUMPDF = m_SUMPDF + m_rightIndexValue;
                    %%% we do not go on the left
                    if  ~m_reachedLeftBorderflag
                        m_leftIndex = m_leftIndex + 1;
                    end;
                end;                
            end;
            taskComplexity = 0;  %%% calculate task complexity as math expectation of PDF
            for i = m_leftIndex : m_rightIndex
                taskComplexity = taskComplexity + m_histXInterval(i) * m_pdfCML(i);
            end;
            taskComplexity = taskComplexity / m_SUMPDF;
            clear('m_SUMPDF','m_leftIndexValue','m_rightIndexValue','m_reachedLeftBorderflag','m_reachedRightBorderflag','i');            
            %%% BUILDIND THRESHOLD SUB ARRAY OF COMPLEXITY INTERVAL
            j = 1;
            m_ThresholdSubList(j) = m_genCML_List( max( find( m_histXInterval(m_leftIndex) >= m_genCML_List ) ) ); 
            for i = max( find( m_histXInterval(m_leftIndex) >= m_genCML_List ) ) : ( max( find( m_histXInterval(m_rightIndex) >= m_genCML_List) ) - 1)
                if ( m_genCML_List(i+1)-m_ThresholdSubList(j) ) > M_CONST_SUBTHRESHOLD_STEP
                    j = j+1;
                    m_ThresholdSubList(j) = m_genCML_List(i+1);
                end;
            end;
            if ( ( max(m_ThresholdSubList) - min(m_ThresholdSubList) ) ~= 0 )
                coefThresholdCompresion = length_ComplexityInterval/( max(m_ThresholdSubList) - min(m_ThresholdSubList) );
            else
                fprintf('m_TDTS.m :: Warning(1) :: coefThresholdCompresion has to be equal +externity, we define it equals 0\n');
                coefThresholdCompresion = 0;
            end;
            clear('m_leftIndex','m_rightIndex','i','j','m_genCML_List');
            m_left_optIndex = 1;
            m_right_optIndex = size(m_ThresholdSubList,2);
            m_optimumList = cell(1,m_right_optIndex); 
            %%%% SAVING INTIAL PARAMETERS now IT IS READY FOR SEARCH FOR OPTIMAL
            save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion');
        end;        
        if isempty(m_optimumList{m_left_optIndex})
            fprintf('m_TDTS.m :: Number of Thresholds to self-tune :: %d\n',m_right_optIndex);
            fprintf('m_TDTS.m :: Self-tuning :: Iteration No. 1\n');
            m_Complexity_PARAMS.Threshold = m_ThresholdSubList(m_left_optIndex);
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                clear('m_results', 'm_results_stat');
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;
            [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, COUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_left_optIndex, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);            
            m_optimumList{m_left_optIndex} = m_results_stat{1,m_left_optIndex}.optimum;
            save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion');
            save(strcat(m_str_outputFileName,'.mat'),'m_results', 'm_results_stat','POUT','COUT','WOUT','TREE','MAXTREELEVEL');    
        end;
        if isempty(m_optimumList{m_right_optIndex})
            fprintf('m_TDTS.m :: Self-tuning :: Iteration No. 2\n');
            m_Complexity_PARAMS.Threshold = m_ThresholdSubList(m_right_optIndex);
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                clear('m_results', 'm_results_stat');
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;            
            [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, COUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_right_optIndex, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);
            m_optimumList{m_right_optIndex} = m_results_stat{1,m_right_optIndex}.optimum;
            save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion');
            save(strcat(m_str_outputFileName,'.mat'),'m_results', 'm_results_stat','POUT','COUT','WOUT','TREE','MAXTREELEVEL');    
        end;
        m_newIndex = floor( (m_right_optIndex-m_left_optIndex)*M_SP_COEF ) + m_left_optIndex;
        if isempty(m_optimumList{m_newIndex})
            fprintf('m_TDTS.m :: Self-tuning :: Iteration No. 3\n* searching for optimum :: [%1.4f;%1.4f]:[%d;%d]\n',m_ThresholdSubList(m_left_optIndex),m_ThresholdSubList(m_right_optIndex),m_left_optIndex,m_right_optIndex);
            m_Complexity_PARAMS.Threshold = m_ThresholdSubList(m_newIndex);
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                clear('m_results', 'm_results_stat');
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;            
            [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, COUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_newIndex, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);
            m_optimumList{m_newIndex} = m_results_stat{1,m_newIndex}.optimum;
            save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion');
            save(strcat(m_str_outputFileName,'.mat'),'m_results', 'm_results_stat','POUT','COUT','WOUT','TREE','MAXTREELEVEL');    
        end;        
        m_deltaTHRIndexes = m_right_optIndex - m_left_optIndex;
        %%%% FIND OPTIMUM IN THE CYCLE
        while (m_deltaTHRIndexes >= 3) %%% while there is more then one element between two, then search for optimal
            if ( m_optimumList{m_newIndex} < min( m_optimumList{m_right_optIndex},m_optimumList{m_left_optIndex} ) )
                if ( m_optimumList{m_left_optIndex} <= m_optimumList{m_right_optIndex} )
                    m_right_optIndex = m_newIndex;   
                else
                    m_left_optIndex = m_newIndex;
                end;
            else
                if ( min( m_optimumList{m_newIndex},m_optimumList{m_left_optIndex} ) <= min( m_optimumList{m_newIndex},m_optimumList{m_right_optIndex} ) )
                    m_right_optIndex = m_newIndex;
                else
                    m_left_optIndex = m_newIndex;            
                end;
            end;
            fprintf('m_TDTS.m :: Self-tuning ::\n* searching for optimum :: [%1.4f;%1.4f]:[%d;%d]\n',m_ThresholdSubList(m_left_optIndex),m_ThresholdSubList(m_right_optIndex),m_left_optIndex,m_right_optIndex);
            m_newIndex = floor( (m_right_optIndex-m_left_optIndex)*M_SP_COEF ) + m_left_optIndex;
            m_Complexity_PARAMS.Threshold = m_ThresholdSubList(m_newIndex);
            if ( exist ( strcat(m_str_outputFileName,'.mat') ) == 2 )
                clear('m_results', 'm_results_stat');
                load(strcat(m_str_outputFileName,'.mat'), 'm_results', 'm_results_stat');
            end;            
            [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, COUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_newIndex, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);
            m_optimumList{m_newIndex} = m_results_stat{1,m_newIndex}.optimum;            
            m_deltaTHRIndexes = m_right_optIndex - m_left_optIndex;
            save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion');
            save(strcat(m_str_outputFileName,'.mat'),'m_results', 'm_results_stat','POUT','COUT','WOUT','TREE','MAXTREELEVEL');    
        end;        
                
        %%% ASSIGN OPTIMUM (SEARCHING FOR MINIMUM )and FIND ITS' INDEX FOR THRESHOLD
        index_opt = 1;
        m_optimum = m_optimumList{index_opt};
        for i = 2 : size(m_optimumList,2)
            if ~isempty(m_optimumList{i})
                if m_optimumList{i} < m_optimum
                    m_optimum = m_optimumList{i};
                    index_opt = i; %%% id wrong the output file is not solid, but consists a different parts
                end;
            end;
        end;
        
        %%%% PRINT RESULTS OF SELF-OPTIMIZATION
        fprintf('m_TDTS.m :: Self-tuning\nOptimal threshold       :  %1.4f\nTask complexity ration  :  %1.4f\nComplexity interval     : [%1.4f;%1.4f]\nSub-complexity interval : [%1.4f;%1.4f]\nCoefficient of compression : %1.4f\n',m_ThresholdSubList(index_opt),taskComplexity,min_cmp_Threshold,max_cmp_Threshold,m_ThresholdSubList(1),m_ThresholdSubList(size(m_ThresholdSubList,2)),coefThresholdCompresion);
        save(m_str_optThresholdFN,'m_avgNumSimCluster','m_stdNumSimCluster','m_avgTreeLevel','m_stdTreeLevel','min_cmp_Threshold','max_cmp_Threshold','m_histXInterval','m_pdfCML','taskComplexity','m_ThresholdSubList','m_left_optIndex','m_right_optIndex','m_optimumList','length_ComplexityInterval','coefThresholdCompresion','m_optimum','index_opt');

        %%% PRINTING OUT FINAL RESULT Of OPTIMUM, code has not been checked before          
        if ~isempty(m_results_stat{1,index_opt})
                fprintf('\nOPTIMAL RESULT / Optimal threshold index : %d\n',index_opt);        
                fprintf('-------------------------------------------------------------\n');
                fprintf('* Average TDTS executing time : %4.4f +/- %4.4f\n* Average Number of Prototypes: %4.4f +/- %4.4f\n* Average Learning rate       : %4.4f +/- %4.4f\n* Average Generalization rate : %4.4f +/- %4.4f\n',m_results_stat{1,index_opt}.avr_ET,m_results_stat{1,index_opt}.std_ET,m_results_stat{1,index_opt}.avr_NPs,m_results_stat{1,index_opt}.std_NPs,m_results_stat{1,index_opt}.avr_LR,m_results_stat{1,index_opt}.std_LR,m_results_stat{1,index_opt}.avr_GR,m_results_stat{1,index_opt}.std_GR);
                fprintf('-------------------------------------------------------------\n');
        else
            fprintf('m_TDTS.m :: Warning(2) :: has not found m_results_stat of optimal index, please check your code\n');
            return;
        end;
        
    otherwise
        error('m_TDTS.m :: Unknown T-DTS mode');
        
end;  % FOR SWITCH OPERATOR :: The end of self-processing mode of TDTS

print_Results_scrpt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SUBFUNCTION LISTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m_results, m_results_stat, TREE, MAXTREELEVEL, POUT, TOUT, WOUT] = m_itT_DTS(m_results, m_results_stat, m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM, dname, m_LearningDBCreating, m_randomizeIncomingData_flag, m_PCA, m_str_Type_ofNormalizing, m_complexityESTList, ind_selected_Complexity, m_thrsh_index, m_iterationNumber, m_Display_PARAMS, m_Complexity_PARAMS, zeroMeaningValue, m_RBFParams4ECTD, resolution_B, m_str_Type_ofDistance, m_DU, m_PU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provide N itteration of T-DTS for fixed range of parameters: Threshold and Complexity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Running N-number of T-DTS Iteration for fixed all parameters
for m_tmp_i = 1 : m_iterationNumber                

    [inPL, inCL, inPG, inCG] = m_uploadDBs(dname, m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
    %%% RERandomizing incoming data DISABLED because it in randome way splited in, so as the all database upload as solid, right now we do not a practical need in this mixing part of code
    %%% This code is not delete because in case of using already splitted on random and generalization part we need random mixing of data
    %if m_randomizeIncomingData_flag
    %    [PL, CL] = m_Dissarange(PL,CL);
    %    [PG, CG] = m_Dissarange(PG,CG);
    %    if m_Display_PARAMS.show_flag
    %        fprintf('m_T_DTS_LrGn.m :: Randomization of Learning and Generalization databases done\n');
    %    end;
    %end;
    fprintf('m_itT_DTS.m :: T-DTS iteration No. %d for %2.4f of %s complexity\n',m_tmp_i, m_Complexity_PARAMS.Threshold, m_complexityESTList{ind_selected_Complexity});
    if m_Display_PARAMS.show_flag
        fprintf('m_T_DTS.m :: m_itT_DTS.m :: Database features ::\n* File-name : %s\n* Learning DB size : %d\n* Generalization DB size : %d\n* Number of vector-attributes : %d\n',dname,size(inCL,2),size(inCG,2),size(inPL,1));
    end;
    
    %%% Now Divide to Simplify
    clear('STR','str','WOUT','POUT','TOUT','m_CML_List');
    STR = [];       % complex struction that accumulates information about decomposition
    m_CML_List = []; % list of the complexity thresholds appeared during decomposition 
    if strcmp(m_DU.method,'none')  %%%% No decomposition go out of procedure and run PU then
        fprintf('m_T_DTS.m :: m_itT_DTS.m :: Note(1) :: No decomposition, as requested\n');
        POUT{1} = inPL;
        TOUT{1} = inCL;
        WOUT = mean(inPL');
        WOUT = WOUT';
        MAXTREELEVEL = 0;
        TREE = [0];
        clear('m_CML_List','str'); 
        m_timeCounter.redec_time = 0;  %%% no re-decomposition , so no time than
    else     
        str='STR'; % additional helpful string that accumulates expression of decomposition
        WOUT = []; % matrix of coordinates of the prototypes after decomposition
        POUT = []; % structure of the decompossed data vectors
        TOUT = []; % structure of the decompossed class labels
        
        tic; % calculating time spent on tree of decomposition buildin 
        [WOUT, POUT, TOUT, TREE, MAXTREELEVEL, m_CML_List] = m_recDBSplit(WOUT,POUT,TOUT,0,inPL,m_getDecMethod_params(m_DU.method),inCL,1,STR,str,m_str_Type_ofDistance,m_Complexity_PARAMS,m_CML_List,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);         
        m_timeCounter.tree_constr_time = toc;           

        clear('m_CML_List','STR','str');
        WOUT = WOUT'; % converting matrix of prototype to standard (like general database of vectors) form  

        if ( strcmp(m_DU.Dec_Param,'Min_dist_to_prototype') )    
            %%%% re-decompose learning database up to not following of the "tree like recursive" paradigm decomposition
            tic; % calculating time spent on re-decomposition
            [POUT, TOUT] = m_splitting(inPL, WOUT, inCL, m_str_Type_ofDistance);
            m_timeCounter.redec_time = toc;   
        else
            m_timeCounter.redec_time = 0;  %%% no re-decomposition
        end;
    end;

    % learning and build structure of ANNs for database clusters and train each of them
    tic;
    for i = 1 : size(TOUT,2) %%% Applying processing units for decompossed database	
        NN{i}.DB_Size = size(TOUT{i},2);
        if ( NN{i}.DB_Size == 0 ) % cought a case when sub-database is empty 
            NN{i}.PU = 'empty_db';
        else    
            NN{i}.PU = m_PU.method;
            NN{i}.net = m_build_sANN(POUT{i},TOUT{i},m_PU.method,m_Display_PARAMS,i,max(inCL)-min(inCL)+1);
        end;
    end;
    m_timeCounter.sANN_building_time = toc; % pick-up the time of learning
    
    % learning of decomposed structure, right now function is called for Learning phase. POUT is decomposed Learning database
    if m_Display_PARAMS.show_flag
        fprintf('m_TDTS.m :: m_itT_DTS.m :: Launching testing for Learning database, then check quality\n');
    end;

    [P, TL] = m_run_sANN(POUT, NN, M_CONST_CLASS_LABEL_UNDETERM);

    % calculate parametes of the quality of classification. It's applied for Learning DB (learning phase)
    tic;
    [m_unsCRate,m_sCRate] = m_QClassification(P,TL,inPL,inCL);
    m_timeCounter.CLearning_time = toc;
    clear('P');
    m_QRates.Learning.Success = m_sCRate * m_CONST_MAX_PRC;
    m_QRates.Learning.Unsuccess = m_unsCRate * m_CONST_MAX_PRC;
    m_QRates.Learning.Undetermined = (1-m_sCRate-m_unsCRate) * m_CONST_MAX_PRC;
    clear('m_unsCRate','m_sCRate');

    % decompose in the same way generalization database
    if ~strcmp(m_DU.method,'none') % catch the case when no decomposition requested, if not of not, then do decomposition
     
        clear('POUTG');
        switch m_DU.Dec_Param
    
            case 'Min_dist_to_prototype'
                tic;
                %%%% CALLING OF RECURSIVE FUNCTION ::: REALIZATION SEE BELOW // CODE CANNOT BE SIMPLYFYID HERE
                [POUTG, TOUTG] = m_splitting(inPG, WOUT, inCG, m_str_Type_ofDistance);
                m_temp_time = toc; 
                clear('TOUTG');
                %%% total time that spent for re-decomposition based on matrix of prototypes, it accumulates for Learning and Generalization database decomposition    
                m_timeCounter.redec_time = m_timeCounter.redec_time + m_temp_time; 
    
            case 'Tree_following'
                POUTG = [];  % structure of the decompossed data vectors
                tic;
                %%%% CALLING OF RECURSIVE FUNCTION ::: REALIZATION SEE BELOW // CODE CANNOT BE SIMPLYFYID HERE              
                POUTG = m_recDBSplitFollowing_tree(POUTG, inPG, inCG, TREE, size(TREE.p,1), m_str_Type_ofDistance);
                m_temp_time = toc;
                m_timeCounter.tree_constr_time = m_timeCounter.tree_constr_time + m_temp_time; %%% total time spent for tree building for Learning and Generalization database
        end;
        clear('m_temp_time');
    else
        %%% no decomposition has been requested, so no decomposition we will do
        POUTG{1} = inPG;
        m_timeCounter.tree_constr_time = 0;    
    end; 

    %%% Learning of decomposed structure. POUTG is decomposed for generalization database
    if m_Display_PARAMS.show_flag
        fprintf('m_itT_DTS.m :: Launching testing for Generalization database, then check quality\n');
    end;
    [P, TG] = m_run_sANN(POUTG, NN , M_CONST_CLASS_LABEL_UNDETERM);

    %%% calculate parametes of the quality of classification. It's applied for generalization phase
    tic;
    [m_unsCRate,m_sCRate] = m_QClassification(P,TG,inPG,inCG);
    m_timeCounter.CGeneralizing_time = toc;
    clear('P');
    m_QRates.Generalizing.Success = m_sCRate * m_CONST_MAX_PRC;
    m_QRates.Generalizing.Unsuccess = m_unsCRate * m_CONST_MAX_PRC;
    m_QRates.Generalizing.Undetermined = (1-m_sCRate-m_unsCRate) * m_CONST_MAX_PRC;
    clear('m_unsCRate','m_sCRate');
    
    %%%% Saving results of each iteration to structure of results                
    m_results{ind_selected_Complexity,m_thrsh_index,m_tmp_i}.Rates = m_QRates;
    m_results{ind_selected_Complexity,m_thrsh_index,m_tmp_i}.NPs = size(WOUT,2);
    m_results{ind_selected_Complexity,m_thrsh_index,m_tmp_i}.Threshold = m_Complexity_PARAMS.Threshold;
    m_results{ind_selected_Complexity,m_thrsh_index,m_tmp_i}.Complexity = m_complexityESTList{ind_selected_Complexity};
    m_results{ind_selected_Complexity,m_thrsh_index,m_tmp_i}.Time_total = m_timeCounter.tree_constr_time +m_timeCounter.redec_time +m_timeCounter.sANN_building_time;
    clear('m_timeCounter');
    
end; %%% end of iteration cycle of T-DTS 
clear('m_tmp_i','m_DU','m_PU','inPL','inPG','TL','TG','inCL','inCG','m_str_Type_ofDistance','m_Complexity_PARAMS','zeroMeaningValue','m_RBFParams4ECTD','resolution_B','m_Display_PARAMS','dname','m_randomizeIncomingData_flag','m_PCA','m_str_Type_ofNormalizing');

%%% calculate statistic
for k_it = 1 : size(m_results,3)
    m_LRa(k_it) = m_results{ind_selected_Complexity,m_thrsh_index,k_it}.Rates.Learning.Success;
    m_GRa(k_it) = m_results{ind_selected_Complexity,m_thrsh_index,k_it}.Rates.Generalizing.Success;
    m_NPs(k_it) = m_results{ind_selected_Complexity,m_thrsh_index,k_it}.NPs;
    m_TT(k_it)  = m_results{ind_selected_Complexity,m_thrsh_index,k_it}.Time_total;
end;
avg_Lr =  mean(m_LRa); avg_Gr =  mean(m_GRa); avg_NPs =  mean(m_NPs); avg_ET =  mean(m_TT);
std_Lr = std(m_LRa)/2; std_Gr = std(m_GRa)/2; std_NPs = std(m_NPs)/2; std_ET = std(m_TT)/2;

m_results_stat{ind_selected_Complexity,m_thrsh_index}.avr_LR  = avg_Lr;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.std_LR  = std_Lr;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.avr_GR  = avg_Gr;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.std_GR  = std_Gr;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.avr_NPs = avg_NPs;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.std_NPs = std_NPs;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.avr_ET  = avg_ET;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.std_ET  = std_ET;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.threshold = m_results{ind_selected_Complexity,m_thrsh_index,k_it}.Threshold;
m_results_stat{ind_selected_Complexity,m_thrsh_index}.optimum = m_optimizationFunction(m_CONST_MAX_PRC, avg_Gr,std_Gr,avg_Lr,std_Lr,avg_ET,std_ET, avg_NPs, std_NPs,m_LearningDBCreating);

clear('m_LearningDBCreating','avg_Lr','avg_Gr','avg_NPs','avg_ET','std_ET','std_Lr','std_Gr','std_NPs','k_it');
clear('m_LRa','m_GRa','m_NPs','m_TT');

function [WOUTnew, POUTnew, TOUTnew, STR, MAXLEVEL,m_CML_List] = m_recDBSplit(WOUTnew, POUTnew, TOUTnew, MAXLEVEL, P,method,C,level,STR,str,m_str_Type_ofDistance,m_Complexity_PARAMS,m_CML_List,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);
% -----------------------------------------------------------------------------------------
% This function is a main recurrent engine of the T-DTS approach. it performs decomposition
% of the database recurrently calling itself when needed, it takes as inputs: vectors P, 
% classes C, decomposition method 'method' the rest of inputs are required for recurrention
% m_indexDistanceType and further the parameters required for subfunctions
% It outputs matrix of final prototypes WOUT and decomposition structure STR.
% STR -structure is STR.n{1,2, M}.p or STR.n{1,2, M}.n{}.n{} : p is an processing unit crd
% STR  - containes decomposiion trace, dadicates to visualization
% WOUTnew - Matrix of Prototypes coordinates, must be transp. before use
% POUTnew{i} - Decomposed vectors' database : i index of dexomposition
% TOUTnew{i} - consists their class labels
% -----------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( (level > MAXLEVEL) | (MAXLEVEL == []) )    
    MAXLEVEL = level;
end;  %%% intializing of the MAXLEVEL // What is the maximal level of tree?

if ( size(C,2) > method.NoN ) %%% split into sub clusters
    W = m_getWPrototype(method,P,C,m_Display_PARAMS);
    [POUT, TOUT] = m_splitting(P, W', C, m_str_Type_ofDistance);
else
    %%% Realization of the special case when requested full decomposition up to the case_level when each sub-cluster has a size less then Number of Prototypes
    if m_Display_PARAMS.show_flag
        fprintf('m_recDBSplit.m :: Note(1) :: Database has reached minimal cluster number\n');
    end;
    m_CML_List(size(m_CML_List,2)+1) = 1; %%% add complexity ratio to the end of list , problem defined as non(1)-complex
    %%% That wierd case has been reached, now stop flag for tree-node, and assign data for node
    eval( strcat(str,' = 0;') );
    WOUTnew(size(WOUTnew,1)+1,:) = mean(P');
    POUTnew{size(POUTnew,2)+1}   = P;
    TOUTnew{size(TOUTnew,2)+1}   = C;   
    return;
end;

%%% Building TREE strucure STR
eval( strcat(str,'.p = W;') );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN RECURSIVE CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : size(TOUT,2)
    if ( size(TOUT{i},2) >= 1 ) %%% check is there is something to think about decomposition        
        %%%% Computing of complexity for new sub-database and taking decide
        [m_flag_TDof_Dec,cm] = m_EC_TD(POUT{i},TOUT{i},method.NoN,m_Complexity_PARAMS,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);
        m_CML_List(size(m_CML_List,2)+1) = cm; %%% add complexity ratio to the end of list
        clear('cm'); %%% clear ration of complexity because we will not use it anymore
        if m_flag_TDof_Dec            
            strnew = strcat(str,'.nnf{',int2str(i),'}');          
			[WOUTnew, POUTnew, TOUTnew, STR, MAXLEVEL, m_CML_List]  = m_recDBSplit(WOUTnew, POUTnew, TOUTnew, MAXLEVEL, POUT{i}, m_getDecMethod_params(method.name), TOUT{i},level+1,STR,strnew,m_str_Type_ofDistance,m_Complexity_PARAMS,m_CML_List,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);   
        else
            strnew_stop = strcat(str,'.nnf{',int2str(i),'}');
            eval( strcat(strnew_stop,' = 0;') );
            WOUTnew(size(WOUTnew,1)+1,:) = W(i,:);  %%% adding at the end of list next coordinate set for nex prototype
            POUTnew{size(POUTnew,2)+1}   = POUT{i}; %%% adding at the end of list next decomposed vector
            TOUTnew{size(TOUTnew,2)+1}   = TOUT{i}; %%% then class label for it
        end;
    else
        strnew_empt = strcat(str,'.nnf{',int2str(i),'}');
        eval( strcat(strnew_empt,' = -1;') );
        if m_Display_PARAMS.show_flag
            fprintf('m_recDBSplit.m :: Note(1) :: Sub-database POUT{i} is empty. No need to call recursive procedure m_recDBSplit\n');
        end;
    end;
end;

function POUT = m_recDBSplitFollowing_tree(POUT, P, T, STR, m_NumProt_onLev, m_str_Type_ofDistance);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This procedure decompose incomming database into sub-database in the same mone how the tree of decomposition %
% has been created on the main recursive complexity estimating based phase of database decomposition           %
% NOTE: IMPLEMENTED ONLY FOR THE TREE WHICH HAS THE SAME NUMBERS OF LEAVES/NODES ON EACH LEVEL               %
% INPUT:                                                                                                       %
% - POUT  input / output decomposed database, creates in dynamic way                                    %
% - P, T - incomind solid database                                                                             %
% - STR - Tree structure                                                                                       %
% - m_NumProt_onLev - number of nodes on each tree-level 2, 3, ...                                             %
% - m_indexDistanceType - parameter that controls which type of distance for decomposition will be applied     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( size(T,2) > m_NumProt_onLev ) %%% When the sub-cluster has size bigger than a matrix of prototypes
    [sPOUT, sTOUT] = m_splitting(P, STR.p', T, m_str_Type_ofDistance);
    for i = 1 : size(sTOUT,2)
        if isstruct( STR.nnf{i} )
           POUT = m_recDBSplitFollowing_tree(POUT, sPOUT{i}, sTOUT{i}, STR.nnf{i}, m_NumProt_onLev , m_str_Type_ofDistance);        
        else
            if ( STR.nnf{i} == 0 )                   % when  it equals -1 it means that this is the case when cluster is/was empty
                if ( size(sTOUT{i},2) >= 1 )         % no need to assign empty clusters
                    POUT{size(POUT,2)+1} = sPOUT{i}; % everything adding at the end of list next decomposed vector                                       
                end;
            end;
        end;
    end;    
else    %%% There is a case when number of vectors in this cluster is less then Number of Prototypes, or probably it's empty
    if (size(T,2) >= 1) %%% When cluster is not empty        
       POUT{size(POUT,2)+1} = P;
    end;    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN SUB FUNCTIONS LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-SUB FUNCTION LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN THIS PART AVALIABLE FUNCTION, that USED BY THE MAIN SUB FUNCTION (function of first call) 

function [m_wrongRate,m_matchRate] = m_QClassification(P,T,m_stdP,m_stdT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function rate of success classification based on standart  
% of DB of class-label and incomping (obtained by some learning
% or generalization process )
% INPUT:
% - P, T - incoming / obtained DB of of vectors and their labels
% - m_stdP, m_stdT - standart
% NOTE P and m_stdP are disordered!
% OUTPUT:
% - Rate of successful classification (gen or learn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrongRate = 0;
matchRate  = 0;
%%% double filtering
if ( size(T,2) ~= size(m_stdT,2) )
    if  ( size(T,2) < size(m_stdT,2) )
        fprintf('m_TDTS.m :: m_QClassification :: Warning(1) :: Discepancy of the size %d\n',(size(m_stdT,2)-size(T,2)) );
    else
        error('m_QClassification :: Error(1) :: output DB is bigger than input');
    end;    
end;

for i = 1 : size(P,2)
    for j = 1 : size(m_stdP,2)
        if ( P(:,i) == m_stdP(:,j) )
            if ( T(i) == m_stdT(j) )
                matchRate = matchRate + 1;
            else
                wrongRate = wrongRate + 1;
            end;
            break;
        end;        
    end;    
end;
m_wrongRate = wrongRate/size(m_stdP,2);
m_matchRate = matchRate/size(m_stdP,2);

function [POUT, TOUT] = m_splitting(X, W, T, m_str_Type_ofDistance);
% ----------------------------------------------------------------------------------------------%
% This function divides the database X and T class labels into subdatabases using vectors of    %
% prototype W (centers) based on approach to pu vector to the cluster where the distance to the % 
% centre is minimal                                                                             %
% INCOMING PARAMETERS:                                                                          %
%  - X - database (matrix of vectors)                                                           %
%  - T - their classes                                                                          %
%  - W - protopyse centers                                                                      %
%  - m_str_Type_ofDistance - selecting type of distance: Eucleadian, Mahalanobis and etc ...    %
% OUTPUT:                                                                                       %
%  - POUT is the structure of splitted databases                                                %
%  - TOUT is the structure of class labels                                                      % 
% ----------------------------------------------------------------------------------------------%

clear('m_distances_Matrix','POUT','TOUT'); 

%%% Filtering incoming data
if (size(W,2) <= 0)
    error('m_splitting.m :: Error(1) :: Database of prototypes W is empty\n');
    return;
end;

if ( (size(W,2) == 1) | (size(X,2) <= 0) )
    POUT = X;
    TOUT = T;
    fprintf('m_splitting.m :: Warnning(1) :: Nothing to decompose\n'); 
    return;
end;

%%% selecting distance type and build matrix of distanceses between vectors and center of clusters
switch m_str_Type_ofDistance
    
  case 'Euclidean' 
     m_distances_Matrix = dist(W',X); 
     
  case 'Manhattan'
     m_distances_Matrix = mandist(W',X);
     
  case 'Mahalanobis'
     error('Mahalanobis distance not realized right now\n');
    
  otherwise
      error('m_splitting.m :: Error(3) :: Unknown distance function, check m_str_Type_ofDistance\n');
end;
clear('m_str_Type_ofDistance');

POUT = [];% structure of the decompossed data vectors
TOUT = [];% structure of the decompossed class labels

[MIND PW] = min(m_distances_Matrix);
clear('MIND','m_distances_Matrix'); 

for i=1 : size(W,2);  %%% splitting process
    indexes = find(PW == i);    %% find indexes of vectors which belong to i prototype ('closer to W(i,:) prototype')
    if ( size(indexes,2) >= 1 ) %% This indexes should exist
        POUT{i} =  X(:,indexes);
	    TOUT{i} =  T(indexes);       
    end;	 
end;     