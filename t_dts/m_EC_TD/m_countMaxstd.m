function m_ratio = m_countMaxstd(X,T);
% ----------------------------------------------------------------------------------------------
% This function finds maximal standard deviation of the SPACE that belongs to diferent classes
% Than it takes max of max and later interpretates (norm) as complexity estimator
% Inputs:
%  -  X, double array of the classification
%  -  T, an array of classes
% Output:
%  -  m_ration
% ----------------------------------------------------------------------------------------------

clear('m_multiClassFDRVector','m_indexFDR','m_classIndex');
for m_classIndex = min(T):max(T)
       indexes = find(T == m_classIndex);;
       m_sumSTD = 0;
       for m_arrayDimension = 1 : size(X,1)
            m_sumSTD = m_sumSTD + abs(std(X(m_arrayDimension,indexes)));               
       end;  %%% end of for m_arrayDimension construction             
       if m_sumSTD < 0
           error('m_countMaxstd.m :: Error(1) :: Coeficient is negative, cannot be used');
       else
            m_multiClassFDRVector(m_classIndex) = m_sumSTD;          
       end;
end;  %%% m_classIndex
m_ratio = 1/(1+max(m_multiClassFDRVector)); %%% Choosing the most complex indicator over multiclassComplexity*components of X(_classIndex,...)