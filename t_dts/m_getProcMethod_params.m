function m_detalized_PU_method = m_getProcMethod_params(m_str_PU_method,m_sizeT,m_totalCN)
% ----------------------------------------------------------------------------------------------------
% This function select the parameters for chooses PU method that will be applied in unique (only for
% this PU methos) for all slaves
% INPUT:
% OUTPUT, structure that consists details of the ap. PU method usage
% ------------------------------------------------------------------------------------------------

m_detalized_PU_method.NoN = 35; %%% as minimum m_totalCN number of neuron(s) will be assigned to first (and maybe one) layer
m_detalized_PU_method.Epo =                 2000; %%% number of epoches of training
m_detalized_PU_method.NoN2=            m_totalCN; %%% second layer - number of classes

switch m_str_PU_method
    
    case {'LVQ1','LVQ2_1'}
        m_detalized_PU_method.LR = 0.01; % learning rate, by default value 0.01 of Param. No. 4, need to be set, be cause of need of use Param. No. 5       
        switch m_str_PU_method            
            case 'LVQ1'
                m_detalized_PU_method.LF = 'learnlv1'; % learning function, this type required for LVQ2.1
            case 'LVQ2_1'
                m_detalized_PU_method.LF = 'learnlv2'; % learning function, this type required for LVQ2.1  
        end;    
        
    case 'RBF'
        %%% Initializing RBF NET, it is need to be set, because of NoN is the last parameter in call function line        
        m_detalized_PU_method.GOAL   =   0;  %%% by default parameter
        m_detalized_PU_method.SPREAD = 1.0;  %%% by default parameter
        m_detalized_PU_method.DF     =  25;  %%% by default parameter        
            
    case {'MLP_FF_BR','Elman_BNBR'}
        %%% Special MLP_FF_BR method for two classes M.Voiry's SAGEM PROBLEM: numbers of epoch 2000, default value from Voiry, first layer of MLP, 20 NNs hat proposed by Voiry M. for test_benchmark and 35 for real problem
        %%% Second range of parameters also proposed by Voiry and we copy them for Elman_BNBR
        m_detalized_PU_method.TFi =  'tansig'; % Transfer function of ith layer, default = 'tansig' / Voiry left it as by default
        m_detalized_PU_method.BTF = 'trainbr'; % Backpropagation network training function, default = 'traingdx', but Voiry modified it to trainbr
        
    case 'MLP_FF_ID'       
        m_detalized_PU_method.ID1 =         0; % first layer delay
        m_detalized_PU_method.ID2 = m_totalCN; % second layer delay
        
end;