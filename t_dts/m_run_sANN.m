function [P, T] = m_run_sANN(POUT, sANN, M_CONST_CLASS_LABEL_UNDETERM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funtion testing trained net in order %
% to obtaine DB, of class labels.           %
% INPUT                                     %
% - POUT - structure of clusters of vectors %
% - sANN - structure of ANN                 %
% OUTPUT                                    %
% - P - solid database of vectors           %
% - T - their class labels                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS FUNCTION CONSISTS ONE SUB FUNCTION FOR INTERPRATATION OF OUTPUT RESULT OF NNet Function

clear('P','T');
P = []; T = [];
%%% Let's get answer and split it in to solid database
for i = 1 : size(POUT,2)
    %%% Check, we do it only for non-empty subdatabases
    if ( ~strcmp(sANN{i}.PU,'empty_db') & (size(POUT{i},2) >= 1) )
        clear('temp_C');
        if isfield(sANN{i}.net.userdata,'m_cluster_label') %%% marked cluster 
                temp_C(1:size(POUT{i},2)) = sANN{i}.net.userdata.m_cluster_label; 
        else 
            %%% apply / recognized simulation of the sANN for i cluster
            % selecting the certaine type of PU and modify parameters of class label output up to the selected method
            switch sANN{i}.PU %%% let's apply trained NN and it will say us what is the class and where    
                           
                case {'LNM','RBF','GRNN'}
                    temp_C = abs(round(sim(sANN{i}.net,POUT{i})));
        
                case {'MLP_FF_BR','Elman_BNwP','Elman_BNBR'}  %%% those methods consist obvious preprocessing so we need to do appopriate convertion
                    warning off;
                    temp_C = m_Matrix2CLdbConvertation( sim( sANN{i}.net,trastd( POUT{i},sANN{i}.net.userdata.MeanP,sANN{i}.net.userdata.StdP ) ), M_CONST_CLASS_LABEL_UNDETERM );
                    warning on;
                    
                otherwise %%% general case
                    warning off;
                    temp_C = m_Matrix2CLdbConvertation( sim( sANN{i}.net,POUT{i}), M_CONST_CLASS_LABEL_UNDETERM );
                    warning on;
            end;  %%% switch end
        end; %%% end of if
        %%% Then check consistance of the output :: VERY IMPORTANT CHECK BEFORE CONCATENATION
        if ( size(POUT{i},2) ~= size(temp_C,2) )
            if ( size(temp_C,2) < size(POUT{i},2) )
                fprintf('m_run_sANN :: NOTE(0) :: %d prototype(s) of the cluster %d - unrecognized\n',(size(POUT{i},2)-size(temp_C,2)),i);
                for cln = size(temp_C,2)+1 : size(POUT{i},2)
                    temp_C(cln) = M_CONST_CLASS_LABEL_UNDETERM;
                end;
            else
                fprintf('m_run_sANN :: NOTE(1) :: %d prototype(s) of the cluster %d - unrecognized\n',size(POUT{i},2),i);
                temp_C(1:size(POUT{i},2)) = M_CONST_CLASS_LABEL_UNDETERM; %%% Marking it as undetermined
            end;
        end;  
        % recognizing done correctly, now concat DBs
        T = [T , temp_C];
        P = [P , POUT{i}];
        %%%% additional check-filter
        if ( size(T,2) ~= size(P,2) )
            error('m_run_sANN :: Error(1) :: Discrepancies between cluster and size labels size found');
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = m_Matrix2CLdbConvertation(m_Matrix,M_CONST_CLASS_LABEL_UNDETERM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts matrix of class label to vector of classes
% The rule of taking decision is concidered to be the following:
% 1. Minimal negative output sim-value will indicate the class label
% 2. If the ration is the minimal value is positive it will indicate unrecognized class
%----------------------------------------------------------------------------------------
% INPUT : Matrix
% OUTPUT: Class label vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('T','m_num_Und_class','i','j');
T = [];
m_num_Und_class = 0;
for j = 1 : size(m_Matrix,2)
    i = find ( min(m_Matrix(:,j)) == m_Matrix(:,j) );
    if ( size(i,1) > 1 )
       T = [T, M_CONST_CLASS_LABEL_UNDETERM];
       m_num_Und_class = m_num_Und_class + 1;
   else
       T = [T, i];
   end;
end; % end of for of j

if (m_num_Und_class > 0)
    fprintf('m_Matrix2CLdbConvertation :: Warning(1) :: %d prototype(s) unrecognized\n',m_num_Und_class);
end;