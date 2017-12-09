function m_ratio = m_FisherDR(X,T);
% --------------------------------------------------------------
% Fisher Disriminant Ratio /Liniear Disriminant Analysis (LDA) /
% Implement for multiclass case:
% 1. Compute fisher ration for each combination of classes
% 2. Choose the maximum(minimum) of the combinations
% Input:
% - X(classIndex,data_componentIndex)
% - T(classIndex) = Number/label of class
% Output:
% - m_ration
% Copyright by Ivan Budnyk 2008 / the 4st of March
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear m_multiClassFDRVector;
m_indexFDR = 0;
for m_classIndex = min(T) : max(T)-1
	for m_classIndex2 = (m_classIndex+1): max(T)
       indexes1 = find(T == m_classIndex);
       indexes2 = find(T == m_classIndex2);
       if (size(indexes1) ~= 0) & (size(indexes2) ~= 0)  %%% This check allows to have holes between class Indexation, for example T = [1 3 4], this procedure is working when no class #2
               % NOTE: Each component of X consist in vector /NOTE MATRIX/ of the data, if you have a need to use X-multidimensional datacube, X(_classIndex, .... ) must be convert to datatype X(_classIndex,_dataIndex) 
               for m_arrayDimension = 1 : size(X,1),  
                   nu1 = mean(X(m_arrayDimension,indexes1));
                   nu2 = mean(X(m_arrayDimension,indexes2));
			       ro1 = var(X(m_arrayDimension,indexes1));
                   ro2 = var(X(m_arrayDimension,indexes2));
                   %%%%%%%%%%%%%%%%% Calculating complexity estimator %%%%%%%%%%%%%%%%%%
                   mFDR_div = ro1+ro2;                     %%%% calculating divisor for Fisher_Complexity_Estimator
                   m_indexFDR = m_indexFDR + 1;            %%%% increasing index for array which consist multiclasses FDR
                   if (mFDR_div == 0)
                       m_multiClassFDRVector(m_indexFDR) = 0; % unepexcted case when sum of sigmas^2 is equal zero
                       fprintf('m_FisherDR.m :: Warning(1) :: Sum of variances is zero, classification problem of the size %d defined as (0)-complex\n',size(T,2));
                       m_ratio = 0;
		       return;
                   else
                       m_multiClassFDRVector(m_indexFDR) = (nu1-nu2)^2/mFDR_div;
                   end;
               end;  %%% end of for m_arrayDimension construction             
       end; %%% end of if_construction              
	end;  %%% m_classIndex2
end;  %%% m_classIndex
m_ratio = min(m_multiClassFDRVector)/(1+min(m_multiClassFDRVector)); %%% Choosing the most complex indicator over multiclassComplexity*components of X(_classIndex,...)