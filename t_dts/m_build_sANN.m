function net = m_build_sANN(P,T,m_str_PUmethod,m_Display_PARAMS,i,m_classLabelsSize);
% ---------------------------------------------------------------------------------------
% This function creates learning ANN models for database of vectors P and class labels T
% INPUT
%  -  P is a sub-database
%  -  m_classLabelsSize is a total number of classes beginning from 1
%  -  T is a database of class-lables
%  -  m_PU is a structure of PU method : type and number of neurons 
%  -  i cluster index
% OUTPUT
%  -  net - created and trained network for sub-cluster
% ---------------------------------------------------------------------------------------------
% CONSISTS ONE SUB-FUNCTION FOR INTERPRATATION OF INCOMING DATA BEFORRE NN-Traning

clear('net');
m_detalized_PU_method = m_getProcMethod_params(m_str_PUmethod,size(T,2),m_classLabelsSize);

if (size(T,2) > 1 )
    
    switch m_str_PUmethod %%% Processing using selected method        
       
        case {'LVQ1','LVQ2_1'}
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: LVQ :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface              
            net = newlvq(minmax(P),m_detalized_PU_method.NoN,m_distrClasses(T),m_detalized_PU_method.LR,m_detalized_PU_method.LF);
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            if (max(T) ~= min(T) ) %%% the cluster consists as minimum 2 classes
                net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            else %%% cluster that consists one pure class
                net.userdata.m_cluster_label = T(1);        
            end;%%%% Train sANN for target-cluster for LVQ
        
        case 'Elman_BN'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: Elman_BN :: Cluster No.%d / %d prototypes\n',i,size(T,2));             
            end; %%% Decided to not apply graphic interface
            net = newelm(minmax(P),[m_detalized_PU_method.NoN,m_detalized_PU_method.NoN2]);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
        
        case 'MLP_CF_GD'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: MLP_CF_GD :: Cluster No.%d / %d prototypes\n',i,size(T,2));             
            end; %%% Decided to not apply graphic interface
            net = newcf(minmax(P),[m_detalized_PU_method.NoN,m_detalized_PU_method.NoN2]);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));    
        
        case 'LNM'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: LNM :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface
            net = newlind(P,T);           
        
        case 'RBF'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: RBF :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface
            warning off;
            net = newrb(P,T);
            warning on;
        
        case 'MLP_FF_GDM'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: MLP_FF_GDM :: Cluster No.%d / %d prototypes\n',i,size(T,2));             
            end; %%% Decided to not apply graphic interface
            net = newff(minmax(P),[m_detalized_PU_method.NoN, m_detalized_PU_method.NoN2]);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
                
        case 'PNN'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: PNN :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface
            net = newpnn(P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            
        case 'GRNN'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: :: GRNN :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface
            net = newgrnn(P,T);
    
        case 'MLP_FF_BR' % Special MLP addapted method for Voiry's problem resolving / has been tested onle for 2 classes labels, we adapt this problem for more classes
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: MLP_FF_BR :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface  
            [pn,meanp,stdp] = prestd(P); 
            net=newff(minmax(pn),[m_detalized_PU_method.NoN, m_detalized_PU_method.NoN2],{m_detalized_PU_method.TFi,m_detalized_PU_method.TFi},m_detalized_PU_method.BTF);   
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;        
            net = train(net,pn,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            %%%% Parameters that will be reguired to do the same transformation but for the testing case
            net.userdata.MeanP = meanp;
            net.userdata.StdP = stdp;
        
        case 'MLP_FF_ID'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: MLP_FF_ID :: Cluster No.%d / %d prototypes\n',i,size(T,2));             
            end; %%% Decided to not apply graphic interface
            net = newfftd(minmax(P),[m_detalized_PU_method.ID1,m_detalized_PU_method.ID2],[m_detalized_PU_method.NoN,m_detalized_PU_method.NoN2]);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;        
            net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
        
        case 'Perceptron'
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: Perceptron :: Cluster No.%d / %d prototypes\n',i,size(T,2));             
            end; %%% Decided to not apply graphic interface
            net = newp(minmax(P),m_detalized_PU_method.NoN2);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;        
            net = train(net,P,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            
        case 'Elman_BNwP' % Special MLP addapted a la Voiry's method for SAGEM defects autoclassifying
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: Elman_BNwP :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface  
            [pn,meanp,stdp] = prestd(P);
            net = newelm(minmax(pn),[m_detalized_PU_method.NoN,m_detalized_PU_method.NoN2]);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;        
            net = train(net,pn,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            %%%% Parameters that will be reguired to do the same transformation but for the testing case
            net.userdata.MeanP = meanp;
            net.userdata.StdP = stdp;
            
        case 'Elman_BNBR' % Special MLP addapted a la Voiry's method for SAGEM defects autoclassifying
            if m_Display_PARAMS.show_flag
                fprintf('m_build_sANN.m :: Elman_BNBR :: Cluster No.%d / %d prototypes\n',i,size(T,2));
            end; %%% Decided to not apply graphic interface  
            [pn,meanp,stdp] = prestd(P);
            net = newelm(minmax(pn),[m_detalized_PU_method.NoN,m_detalized_PU_method.NoN2],{m_detalized_PU_method.TFi,m_detalized_PU_method.TFi},m_detalized_PU_method.BTF);
            net = init(net);
            net.trainParam.epochs = m_detalized_PU_method.Epo;        
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;        
            net = train(net,pn,m_CLdb2MatrixConvertation(T,m_classLabelsSize));
            %%%% Parameters that will be reguired to do the same transformation but for the testing case
            net.userdata.MeanP = meanp;
            net.userdata.StdP = stdp;

        
        otherwise
            error('m_build_sANN.m :: Error(1) :: No such PU method, please check incoming parameters');
    end; %end of switch
    
else %%% cluster consists one vector of prototype
    net.userdata.m_cluster_label = T(1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SUB FUNCTION LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m_Matrix = m_CLdb2MatrixConvertation(T,m_classLabelsSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function convert class db to Matrix forma
% OUTPUT : Matrix
% INPUT: Class label vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_Matrix = ones(m_classLabelsSize,size(T,2));

for j = 1 : size(T,2)
    for i = 1 : m_classLabelsSize
        if (T(j) == i)  
            m_Matrix(i,j) = -1; %%% class labeling
            break; %%% we have picked-up the case to indicate the class so no need to go to the end of cycle "for j"
        end;
    end; % enf of for of i
end; % end of for of j