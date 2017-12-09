function cml_ratio = m_RBF_ZISC_fusion(X,T,dstn_type,MIF,n_iteration,m_pol_pow,m_Display_PARAMS);
% ---------------------------------------------------------------------------------------------------------
% Simulated ZISC complexity estimator, which is based on the idea of esimating complexity for IBM ZISC-036
% This main approach has been applied using RBF MATLAB function, instead of fixed ZISC RBF net type
% After calculating complexity ration, we normalized it using information about reached performance
% Inputs:
%  -  X, double array of the classification data
%  -  T, an array of the target(numbers) of classes
%  -  dstn_type, type of distance to calculate the difference between two prototypes
%  -  n_iteration, number of iteratin to get average ratio of complexity
%  -  m_pol_pow, polynome power for aproximation of the tail of the basic complexity ratio behaviour
% ---------------------------------------------------------------------------------------------------------
% THIS PROCEDURE CALL MAIN SUB PROCEDURE m_ZISC_CE

if ( MIF <= 0  ) %%% verifying incoming MIF if non accepteable set parameters by default
    m_MIF = m_distance(min(X'),max(X'),dstn_type);
else
    m_MIF = MIF;
end;

for i = 1 : n_iteration
    m_ratio(i) = 1-m_ZISC_CE(X,T,m_MIF,m_pol_pow,dstn_type,m_Display_PARAMS.show_flag);
end;
clear('i','m_MIF');
%% special code for testing this estimator and it's paraleters, required for research and before industreal use
fprintf('M_AVR_RATIO = %6.4f\nM_STD_RATIO=%6.4f\n',mean(m_ratio),std(m_ratio)); 
cml_ratio = mean(m_ratio);
return;

function m_ratio = m_ZISC_CE(P,C,MAXIF,M_CONST_POLPOW,M_STR_DISTANCE,m_show_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This procedure is a simulation of ZISC RBF base algorithm of classification / clustering    %
% Besides clustering, we used this approach for classification complexity estimating          %
% INPUT :                                                                                     %
%   P - database                                                                              %
%   C - class labels                                                                          %
%   MAXIF - maximal influence field                                                           %
%   M_STR_DISTANCE - selection of the type of distance                                        %
%   M_CONST_POLPOW - polynome power which approximates function                               %
% OUPUT :                                                                                     %
%   m_ratio _ ratio of classificational complexity                                            %
% (c) Ivan BUDNYK 2008                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS PROCEDURE CALLS SUB PROCEDURE m_distance to calculate the distance beetween two vectors

clear('m_listIND_USED_data','m_ZISCNNStructure','m_numNeurons');

%%% a counter of neurons in m_ZISCNNStructure
m_numNeurons = 1;                                                     

%%% drop random drop as a centre of the first cluster 
i = round( rand*( size(C,2)-1 ) ) +1;
%%% INIT FIRST NEURON OF ZISCNNStructure
m_ZISCNNStructure(m_numNeurons).CENTRE = P(:,i);
m_ZISCNNStructure(m_numNeurons).CAT    = C(i);
m_ZISCNNStructure(m_numNeurons).ROI    = MAXIF;
m_ZISCNNStructure(m_numNeurons).PROT(:,1)= P(:,i);

%%% mark the data in P as such as used
m_listIND_USED_data(1) = i; 

while (  size(m_listIND_USED_data,2) < size(C,2) ) %%% continue clusterizing if there is still a prototype that not marked as used
    %%% select all prototypes that belongs to MAXIF of the new neuron and that prototypes that have not been checked
    m_index_ofUsedPrototypes = 1;
    %%% calculate the number of clusters that are in ROI
    m_counter_ofINCluster = 0;

    for m_index_ofPrototypes = 1 : size(C,2)
        
       if ( m_index_ofPrototypes ~= m_listIND_USED_data(m_index_ofUsedPrototypes) )
          %%% calculate the distance between centre of neuron and selected protype
          m_getDistance = m_distance(m_ZISCNNStructure(m_numNeurons).CENTRE,P(:,m_index_ofPrototypes),M_STR_DISTANCE);
          if ( m_getDistance <= MAXIF)
              m_counter_ofINCluster = m_counter_ofINCluster + 1;
              m_bigCluster_Prototype(:,m_counter_ofINCluster)= P(:,m_index_ofPrototypes);
              m_bigCluster_Category(m_counter_ofINCluster) = C(m_index_ofPrototypes);
              m_bigCluster_Index(m_counter_ofINCluster)    = m_index_ofPrototypes;
              m_bigCluster_Distances(m_counter_ofINCluster)= m_getDistance;
          end;
       else
           if ( m_index_ofUsedPrototypes < size(m_listIND_USED_data,2) )
                m_index_ofUsedPrototypes = m_index_ofUsedPrototypes + 1;
           end;
       end;
    end;
    clear('m_getDistance','m_index_ofUsedPrototypes','m_index_ofPrototypes');
    
    if (m_counter_ofINCluster ~= 0) %%% if inside the cluster there is a prototype except neuron then narrow ROI of the big cluster
    
        [m_sorted_array_ofDistances, INDEXES] = sort(m_bigCluster_Distances); %%% increase ROI begins from the shortest one
        for m_index_ofPrototypes = 1 : size(m_sorted_array_ofDistances,2)
            %%% when it's still the same category
            if ( m_bigCluster_Category(INDEXES(m_index_ofPrototypes)) == m_ZISCNNStructure(m_numNeurons).CAT )
                %%% than continue build cluster of ZISCNNStructure
                % increase ROI by taking next biger one
                m_ZISCNNStructure(m_numNeurons).ROI = m_sorted_array_ofDistances(m_index_ofPrototypes);
                % add prototype to cluster
                m_ZISCNNStructure(m_numNeurons).PROT(:,m_index_ofPrototypes+1)= m_bigCluster_Prototype(:,INDEXES(m_index_ofPrototypes));
                % mark prototype as used
                m_listIND_USED_data(size(m_listIND_USED_data,2)+1) = m_bigCluster_Index(INDEXES(m_index_ofPrototypes));    
            else
                %%% exit, because we cannot add data to cluster because it's another category
                break;
            end;
        end;    
        m_listIND_USED_data = sort(m_listIND_USED_data);
        clear('m_bigCluster_Prototype','m_bigCluster_Category','m_bigCluster_Index','m_bigCluster_Distances','m_index_ofPrototypes','INDEXES','m_sorted_array_ofDistances');     
        
    end;
    
    clear('m_counter_ofINCluster');
    
    if ( size(m_listIND_USED_data,2) < size(C,2) ) %%% if there is still someting to clusterize let's cluster it, 
                                                   %%% let's add new neuron to m_ZISCNNStructure                                                   
        m_numNeurons = m_numNeurons + 1; %%% add new neuron
        %%% new index of random cluster, but where exatcly it locates in alread clusterized database we have to find
        m_newIndx_randCluster = round( rand*( size(P,2)-size(m_listIND_USED_data,2)-1 ) ) +1; 
        %%% find the centre for new neuron in a terms of index of marked P-database
        k = 1;
        
        for i = 1 :size(C,2)
            
            if (m_listIND_USED_data(k) ~= i )
                m_newIndx_randCluster = m_newIndx_randCluster - 1; %%% here we've catched a prototype not used                
                if (m_newIndx_randCluster == 0) %%% case of exist from the cycle when random number in term of marked data-P vectors belongs to sub-interval of hole size that consist marked vectors
                    break;
                end;    
            else
                k = k + 1; %%% and here we catched marked one, go ahead anyway
                if (size(m_listIND_USED_data,2) < k)
                    i = m_newIndx_randCluster + m_listIND_USED_data(size(m_listIND_USED_data,2));
                    break;
                end;
            end;
        end;
        %end;
        clear('k','m_newIndx_randCluster');
        
        %%% added it to ZISC Neural Network Structure
        m_ZISCNNStructure(m_numNeurons).CENTRE = P(:,i);
        m_ZISCNNStructure(m_numNeurons).CAT    = C(i);
        m_ZISCNNStructure(m_numNeurons).ROI    = MAXIF;
        m_ZISCNNStructure(m_numNeurons).PROT(:,1)= P(:,i);        
        %%% marked it
        m_listIND_USED_data(size(m_listIND_USED_data,2)+1) = i;
        m_listIND_USED_data = sort(m_listIND_USED_data);
        
        clear('i');        
        
    end; %%% end of adding new neuron
    
end; %%% end of clusterizing procedure

clear('m_Qarray','m_sum_ofPrototype');
m_sum_ofPrototype = 0;
for i = 1 : size(m_ZISCNNStructure,2)
    m_sum_ofPrototype = m_sum_ofPrototype + size(m_ZISCNNStructure(i).PROT,2);
    m_Qarray(size(m_ZISCNNStructure,2)-i+1) = i/m_sum_ofPrototype;
end;

clear('m_sum_ofPrototype');

if ~m_show_flag %%% verify do we need print out warning in case of bad approximation
    warning off;
end;
[p,S,mu] = polyfit( (1 : size(m_Qarray,2)), m_Qarray, M_CONST_POLPOW);
if ~m_show_flag %%% related to previouse code
    warning on;
end;

md2Q = polyder(polyder(p));
roots_md2Q = roots(md2Q);

for i = 1 :size(roots_md2Q,1)
    if isreal( roots_md2Q(i) )
        m_ratio = polyval(p,roots_md2Q(i),S,mu);        
        if ( ( m_ratio < 0 ) | ( m_ratio > 1 ) )
            fprintf('m_ZISC_CE.m :: Warning(0) :: Complexity was out of range, ratio has been set as maximal\n');
            m_ratio = m_Qarray(1);             
        end;
        return;
    end;
end;
if m_show_flag
    fprintf('m_RBF_ZISC_fusion.m :: m_ZISC_CE :: Warning(1) :: d2Q(x)=0 has not been found, ratio has been set as maximal\n');
end;
m_ratio = m_Qarray(1);

%%%%%%%%%%%%%%%%%%% SUBFUNCTION OF THE DISTANCE CALCULATION %%%%%%%%%%%%%%%%%%%
function mydistance = m_distance(V1,V2,M_STR_DISTANCE)

clear('mydistance','tmpIndex');

switch M_STR_DISTANCE
    
    case 'L1' %%% City block / Manhat. distance , equvalent of pdist(X,'CityBlock')
        mydistance = 0;
        for tmpIndex = 1 : size(V1,1)
            mydistance = mydistance + abs( V1(tmpIndex)-V2(tmpIndex) ) ;
        end;

    case 'LSUP' %%% Hyper-cube influence field
        mydistance = abs( V1(1)-V2(1) );
        for tmpIndex = 2 : size(V1,1)
            if ( abs(V1(tmpIndex)-V2(tmpIndex)) > mydistance )
                mydistance = abs( V1(tmpIndex)-V2(tmpIndex) );
            end;
        end;
        
    case 'EUCL' %%% Euclidian distance
        mydistance = 0;
        for tmpIndex = 1 : size(V1,1)
            mydistance = mydistance + ( V1(tmpIndex)-V2(tmpIndex) )^2 ;
        end;
        mydistance = sqrt(mydistance);
        
    otherwise   %%% 
        error('m_RBF_ZISC_fusion :: m_distance :: Error(1) :: Type of the distance is not undefine');
end;