function W = m_getWPrototype(method,P,T,m_Display_PARAMS);
% ---------------------------------------------------------------------------
% This function finds the prototypes for given database P and class labels C, 
% using methods, define as input parameter
% INPUT
%  -  P - database
%     method.NoN - number of neurons
% OUTPUT
%  -  W - matrix of the centers of decomposition , must be Transported
% ---------------------------------------------------------------------------

% filtering incomming data
if ( size(P,2) <= 1 )
    error('m_getWPrototype.m :: Error(1) :: No database for decomposing, this case should be cought before NN decomposition usage!');
end;
% 
if ( method.NoN <= 1 )
    fprintf('m_getWPrototype.m :: Warning(1) :: method.NoN should be more than 2, it is reinitialized method.NoN = 2\n');
end;

clear('net','i','m_SOM_grid');
switch method.name %% processing decomposition following selected method type
    
    case 'CNN'
        if m_Display_PARAMS.show_flag
            fprintf('m_getWPrototype :: CNN started\n');
        end;
        net = newc(minmax(P),method.NoN);
        net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
        net = init(net);
        net = train(net,P);
        
    case 'SOM'        
        if m_Display_PARAMS.show_flag
            fprintf('m_getWPrototype :: SOM started\n');
        end;
        net = newsom(minmax(P),method.Grid);
        net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
        net = init(net);
        net = train(net,P);
        
    case {'LVQ1','LVQ2_1'}        
        
        % Created network
        M_CONS_LVQ_LR_DEF = 0.01;
        
        switch method.name
            
            case 'LVQ1'
                if m_Display_PARAMS.show_flag
                    fprintf('m_getWPrototype :: LVQ1 started\n');
                end;
                net = newlvq(minmax(P),method.NoN,m_distrClasses(T),M_CONS_LVQ_LR_DEF,'learnlv1'); 
                
            case 'LVQ2_1'
                if m_Display_PARAMS.show_flag
                    fprintf('m_getWPrototype :: LVQ2_1 started\n');
                end;                
                net = newlvq(minmax(P),method.NoN,m_distrClasses(T),M_CONS_LVQ_LR_DEF,'learnlv2');         
        end;        
        %%% Train it
        if (max(T) ~= min(T) )
            net.trainParam.show = m_Display_PARAMS.show_learningSlaves_batchSize;
            net = init(net);
            net = train(net,P,ind2vec(T)); %%%% Train sANN for target-cluster for LVQ
        else
            fprintf('m_getWPrototype :: Warning(1) :: LVQ detected one class cluster, decomposition skipped\n');
        end;
               
    case 'none'
		% No decomposition, skip doing anything
        
end; %%% switch end

% Get the prototypes, like in example code
W = net.IW{1,1};