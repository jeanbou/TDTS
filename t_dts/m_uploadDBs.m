function [inPL, inCL, inPG, inCG] = m_uploadDBs(dname, m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upload of databases, preprocess them: normalizied, split, provide PCA/T and gives ouput separate learning and generalization parts     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS FUNCTION USE SUBFUNCTION m_randData to randomize incoming data in chaotic way on request

% Upload data from mat file and check basic properties / filtering
    if (exist(strcat(dname,'.mat')) == 2 )
            load(strcat(dname,'.mat'));
    else
            error('m_uploadDBs.m :: Error(1) :: Datafile has not been found');
    end;

    %%% if it is solid database , so split itupload database for learning and generalization
    if ( m_randomizeIncomingData_flag & ( exist('P') & exist('C') ) )
        
        if ~( exist('P') & exist('C') )
            error('m_uploadDBs.m :: Error(2) :: One of the database P or C (or both) is (are) absent');
        end;

        if (size(P,2) ~= size(C,2))
            error('m_uploadDBs.m :: Error(3) :: Discrepancies of size between DB vectors and DB class-labels');
        end;

        if ( min(C) < 1 )
            error('m_uploadDBs.m :: Error(4) :: Class label C should begin from 1, please re-label class vector');
        end;
        
         if ( min(C) == max(C) )
            error('m_uploadDBs.m :: Error(5) :: Nothing to classifying, incoming data belongs to one class');
        end;
        
        for i = min(C) : max(C)
            check_index = find(C == i);
            if ~check_index
                error('m_uploadDBs.m :: Error(6) :: In DB of class labels exists a hole between class numerization, please re-label this DB');
            end;
        end;
    
        %%%% Extracting sub-database from solid learning database following incoming parameters of m_LearningDBCreating structure
        if ( (m_LearningDBCreating.percentage > m_CONST_MAX_PRC) | (round(m_LearningDBCreating.percentage) <= 0) )
            error('m_uploadDBs.m :: Error (7) :: Percentage of learning database size is out of range');
            return;
        end;
        clear('inPL','inCL','inPG','inPL');
        switch m_LearningDBCreating.mode
    
            case 'Random'
                m_vectorOFRandomIndexes = m_randData(round(m_LearningDBCreating.percentage*size(C,2)/m_CONST_MAX_PRC),size(C,2));
                m_vectorOFRandomIndexes = sort(m_vectorOFRandomIndexes);
                jl = 1;
                jg = 1;
                for i = 1 : size(C,2)
                    if ( jl <= size(m_vectorOFRandomIndexes,2) & ( i == m_vectorOFRandomIndexes(jl) ) )
                        inCL(jl) = C(i);
                        inPL(:,jl) = P(:,i);
                        jl = jl+1;
                    else
                        inCG(jg) = C(i);
                        inPG(:,jg) = P(:,i);
                        jg = jg+1;
                    end;
                end;
                clear('jg','jl','m_vectorOFRandomIndexes','i');
        
            case 'Balanced'
                        
                inCL = [];
                inPL = [];
                inCG = [];
                inPG = [];
                for i = min(C):max(C)
                    m_Indexes = find(i == C);            
                    if ( size(m_Indexes,2) == 0 )
                        error('m_uploadDBs.m :: Error (8) :: No random element of one of the smallest classes has not been choosen, please increase (decrease) percentage');
                        return;
                    else
                        m_vectorOFRandomIndexes = m_randData(round(m_LearningDBCreating.percentage*size(m_Indexes,2)/m_CONST_MAX_PRC),size(m_Indexes,2));
                        m_vectorOFRandomIndexes = sort(m_vectorOFRandomIndexes);
                        jl = 1;
                        for i = 1 : size(m_Indexes,2)
                            if ( jl <= size(m_vectorOFRandomIndexes,2) & ( i == m_vectorOFRandomIndexes(jl) ) )
                                inCL(size(inCL,2)+1) = C(m_Indexes(i));
                                inPL(:,size(inPL,2)+1) = P(:,m_Indexes(i));
                                jl = jl+1;
                            else
                                inCG(size(inCG,2)+1) = C(m_Indexes(i));
                                inPG(:,size(inPG,2)+1) = P(:,m_Indexes(i));                
                            end;
                        end;
                    end;
                end;
                clear('m_vectorOFRandomIndexes','jl','m_Indexes','i');
                
            otherwise
                error('m_uploadDBs.m :: Error (9) :: Wrong mode of database splitting');       
        end; %%% end of switch two variants of database splitting
        clear('P','C');
    
    end; % end of the if where P an C splitting into inPL, inPG, inCL and inCG, now else
        
    switch m_str_Type_ofNormalizing
    
        case 'Norm_LinTrans_EC'
            for i = 1 : size(inPL,1)
                inPL(i,:) = (inPL(i,:)-min(inPL(i,:)))/(max(inPL(i,:))-min(inPL(i,:)));
            end;
            for i = 1 : size(inPG,1)
                inPG(i,:) = (inPG(i,:)-min(inPG(i,:)))/(max(inPG(i,:))-min(inPG(i,:)));
            end;
            clear('i');
            
        case 'Norm_Statistical_EC',
            if ~m_PCA.flag
                [inPL,meanp,stdp] = prestd(inPL);
                [inPG,meanp,stdp] = prestd(inPG);
                clear('meanp','stdp');
            else
                fprintf('m_uploadDBs.m :: executing Norm_Statistical_EC normalization mode has been skipped, because PCA/T will be applied with the same normalization\n');
            end;
          
    end; %%% end of switch of normalization DBs techniques
    
    if m_PCA.flag

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Principal Component Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [inPL,m_PCA.meanp,m_PCA.stdp] = prestd(inPL);
        [inPL,m_PCA.transMat] = prepca(inPL,m_PCA.per_elmnt_cmps/m_CONST_MAX_PRC); 
        %%% check how the PCA has been done regardering the database size
        if isempty(inPL)
            error('m_uploadDBs.m :: Error (10) :: PCT removed all components in Learning database. Now it is empty, please reduce PCA rate');
            return;
        end;        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Principal Component Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % After the network has been trained, m_PCA.transMat matrix should be used to transform any future 
        % inputs that are applied to the network. It effectively becomes a part of the network, just like
        % the network weights and biases. 
        % If prepca has been used to preprocess the training set data before, then whenever the trained 
        % network is used with new inputs they should be preprocessed with the transformation matrix that
        % was computed for the training set. This can be accomplished with the routine trapca.
        inPG = trapca(trastd(inPG,m_PCA.meanp,m_PCA.stdp),m_PCA.transMat);
    end;
    
    %%%% filtering of incoming inPL, inPG, inCL and inCG that already splitted on learning and generalization databases or directly downloaded 
    
    if ~( exist('inPL') & exist('inPG') & exist('inCG') & exist('inCL') )
        error('m_uploadDBs.m :: Error(11) :: One or more Databases of inPL, inPG, inCL and inCG is (are) absent');
    end;

    if ( (size(inPL,2) ~= size(inCL,2)) | (size(inPG,2) ~= size(inCG,2)) )
        error('m_uploadDBs.m :: Error(12) :: Discrepancies of size between DBs inPL and inCL or of vectors or inPG and inCG');
    end;

    if ( (size(inPL,1) ~= size(inPG,1)) | (size(inCG,1) ~= size(inCL,1)) )
        error('m_uploadDBs.m :: Error(13) :: Discrepancies of the vector size in DBs or label-class DB');
    end;

    if ( (min(inCL) < 1) | (min(inCG) < 1) )
        error('m_uploadDBs.m :: Error(14) :: Class inCL or inCG label-class DBs should begin from 1, please re-label class vector');
    end;

    if ( (size(inCL,2) == 0) | (size(inCG,2) == 0) )
        error('m_uploadDBs.m :: Error(15) :: Class inCG or inCL DBs class label is (are) empty');
    end;

    if ( (size(inPG,2) == 0) | (size(inPL,2) == 0) )
        error('m_uploadDBs.m :: Error(16) :: DB(s) inPG or inPL is (are) empty');
    end;
    
    for i = min(inCL) : max(inCL)
        check_index = find(inCL == i);
        if ~check_index
            error('m_uploadDBs.m :: Error(17) :: inCL DB of class labels exist a hole between class numerization, please re-label this DB');
        end;
    end;

    for i = min(inCG) : max(inCG)
        check_index = find(inCG == i);
        if ~check_index
            error('m_uploadDBs.m :: Error(18) :: inCG DB of class labels exist a hole between class numerization, please re-label this DB');
        end;
    end;
    clear('i','check_index');
    
    if (  ( max(inCG)-min(inCG) ) ~= ( max(inCL)-min(inCL) ) )
          error('m_uploadDBs.m :: Error (19) :: Outgoing splitted DBs of class labels are disbalanced, please check you learning splitting DB rate or incoming data');
    end; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m_vectorsToZISC = m_randData(MAX_NUM_VEC_FR_LEARNING,m_INTERVAL);
% ------------------------------------------------------------------------------
% This function retunds the array size MAX_NUM_VEC_FR_LEARNING of random indexes
% taken from the interval [1:m_Interval]
% ------------------------------------------------------------------------------

if nargin < 2
  error('m_randData.m :: Error(1) :: Check incoming parameters, please');
  return
end

clear m_tempIntIndexArray;
clear m_vectorsToZISC;

m_tempIntIndexArray = [1 : m_INTERVAL];

for i= 1 : MAX_NUM_VEC_FR_LEARNING
    k = round(rand*(m_INTERVAL-i))+1;
    m_vectorsToZISC(i) = m_tempIntIndexArray(k);
    for j = k : m_INTERVAL - i
        m_tempIntIndexArray(j) = m_tempIntIndexArray(j+1);  
    end;
end;