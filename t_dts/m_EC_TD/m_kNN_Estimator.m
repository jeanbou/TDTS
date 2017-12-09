function m_ratio = m_kNN_Estimator(X,T,k);
% --------------------------------------------------------------
% m_kNN_Estimator estimates Bayes error, building asymptotically curve 
% Implement for multiclass case:
% Input:
% - T(classIndex) = Number/label of class
% - k resolution parameter  
% Output:
% - m_ration
% Copyright by Ivan Budnyk 2008 / the 9th of March
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (k < 2)
    fprintf('m_kNN_Estimator :: Note(1) :: k - resolution is less than 2, k set to equal 2\n');
    k = 2;
end; % check of the incoming parameter B / k

if (size(T,2) < k) %%% It's mandatory check when k will be incoming parameter
    fprintf('m_kNN_Estimator :: Note(2) :: k - resolution is largen than data size, k set to max number of the components\n');
    k = size(T,2);	
end;

m_maxClassNumber = max(T);
m_minClassNumber = min(T);

for i_attribute = 1:size(X,1)
    XmaxAtt(i_attribute) = max(X(i_attribute,:));
    XminAtt(i_attribute) = min(X(i_attribute,:));
end;

clear tmp_array_ofRadiuses;
clear Tnew;
flag_nonZeroRadius = 0;

for i_attribute = 1:size(X,1)
    if ( XmaxAtt(i_attribute)-XminAtt(i_attribute) ~= 0 )
        flag_nonZeroRadius = 1;
        tmp_array_ofRadiuses(i_attribute) = XmaxAtt(i_attribute)-XminAtt(i_attribute);  %%%% we take the minimal radius over all components (directions) of X
    end;
end;

if ~flag_nonZeroRadius
	fprintf('m_kNN_Estimator :: Warning(1) :: components in all direction are concentrated in one dot, problem predefined as 1(non)-complex\n');
    m_ratio = 1;
end;

clear Tnew;
Tnew = T;
if (m_minClassNumber <= 0) %%% find this classes where labels of classes equal zero, and relable it. It is need for correct predefining of the claster
    indexes = find(T == 0);
    m_maxClassNumber = m_maxClassNumber + 1;
    Tnew(indexes) = m_maxClassNumber;
    m_minClassNumber = min(Tnew);
end;

for i_attribute = 1:size(X,1) %%% random choose of the initial centre of the 1st cluster // random choose of the X0 for all components in all directions
    X0(i_attribute) = 0.5*(XmaxAtt(i_attribute)+XminAtt(i_attribute));    
end;
X0 = X0'; %%% There is no need to do such traformation however in variable viewer it will looks like for X

m_deltaR = min(tmp_array_ofRadiuses)/k;
clear tmp_array_ofRadiuses;  %%% no need in this variable anymore
% Calculating maximal possibla radius from XO dot
clear max_dist_by_Coordinate;
for i = 1:size(X,1)
    max_dist_by_Coordinate(i) = max( abs(X(i,:)-X0(i)) );
end;
max_pos_rad = max(max_dist_by_Coordinate);
clear max_dist_by_Coordinate;

clear m_Clusters; %%% Array which consists the cluster
N = 1;
m_tmp_radius = 0;
while ( m_tmp_radius^2 <= (max_pos_rad+m_deltaR)^2 )
    m_componentNumber = 1;
    for j =1 : size(X,2)
        SSx = 0;
        for i = 1 : size(X,1)
            SSx = SSx + (X(i,j)-X0(i))^2;            
        end;
        if ( (SSx < (m_tmp_radius+m_deltaR)^2) & (SSx >= m_tmp_radius^2) )
                m_Clusters(N,m_componentNumber) = Tnew(j);  %%% we put in cluster a number of class of the component 
                m_componentNumber = m_componentNumber + 1;
        end;
    end;
    N = N + 1; %%% new cluster building
    m_tmp_radius = m_tmp_radius + m_deltaR;
end; %%% Main cycle that control cluster creating process

if ~exist('m_Clusters') %%% Partitioning parameter k, and initial random X0 is so that m_Clusters has no been created, so in this case m_Cluster will be solid and will consist all segment
    m_Clusters = Tnew;
end;

sum_4N= 0;
N_Parameter = 0;
for i = 1:size(m_Clusters,1)
    min_classNum_inCluster = min(m_Clusters(i,:));
    max_classNum_inCluster = max(m_Clusters(i,:));
    if ( min_classNum_inCluster~=max_classNum_inCluster )
        max_NumClasses = 0;
        num_Instance = 0;
        for i_classNum = m_minClassNumber:m_maxClassNumber  %%% We take to account only classes in claster than not equal zero. min_classNum_inCluster can consist zero
          if (i_classNum ~= 0)%%% exclude 0-classes, zero in C array mean no X in the subspace
             indexes = find(m_Clusters(i,:)==i_classNum);            
             if ( size(indexes,2) >= max_NumClasses )
                    max_NumClasses = size(indexes,2);
             end;
             num_Instance = num_Instance + size(indexes,2);
          end; % 0-class check if   
        end; % by attributes
        sum_4N = sum_4N + max_NumClasses/num_Instance;
        N_Parameter=N_Parameter+1;
    end;  % if which the case when in choosen cluster existes differents patterns which belongs different patterns!  
end;  % by clusters attributes
m_ratio = sum_4N/N_Parameter; % Choosing the most complex indicator over multiclassComplexity*components of X(_classIndex,...)