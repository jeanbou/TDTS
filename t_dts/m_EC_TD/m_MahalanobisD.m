function m_ratio = m_MahalanobisD(X,T,err_ZeroDet,m_Display_PARAMS);
% --------------------------------------------------------------
% Mahalanobis distance
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
m_maxClassNumber = max(T);
m_minClassNumber = min(T);
clear m_multiClassFDRVector;
m_indexFDR = 0;
for m_classIndex = m_minClassNumber:m_maxClassNumber,
	for m_classIndex2 = (m_classIndex+1):m_maxClassNumber,
       indexes1 = find(T == m_classIndex);
       indexes2 = find(T == m_classIndex2);
       if (size(indexes1) ~= 0) & (size(indexes2) ~= 0)  %%% This check allows to have holes between class Indexation, for example T = [1 3 4], this procedure is working when no class #2
               % NOTE: Each component of X consist in vector /NOTE MATRIX/ of the data, if you have a need to use X-multidimensional datacube, X(_classIndex, .... ) must be convert to datatype X(_classIndex,_dataIndex)
               clear('XC1','XC2','XC1mean','XC2mean');
               XC1 = X(:,indexes1);
               XC2 = X(:,indexes2);               
               %---------- Calculating of the mean vector
               for m_ComponentNumber=1:size(XC1,1)
                   XC1mean(m_ComponentNumber) = mean(XC1(m_ComponentNumber,:));
               end;
               
               for m_ComponentNumber=1:size(XC2,1)
                   XC2mean(m_ComponentNumber) = mean(XC2(m_ComponentNumber,:));
               end;
               
               %   MAHAL Mahalanobis distance.
               %   MAHAL(Y,X) gives the Mahalanobis distance of each point in
               %   Y from the sample in X. 
               %   
               %   The number of columns of Y must equal the number
               %   of columns in X, but the number of rows may differ.
               %   The number of rows must exceed the number of columns
               %   in X.

               %   B.A. Jones 2-04-95
               %   Copyright 1993-2000 The MathWorks, Inc. 
               %   $Revision: 2.8 $  $Date: 2000/05/26 18:53:01 $
               XC1mean = XC1mean';
               XC2mean = XC2mean';
               [rx,cx] = size(XC1mean);
               [ry,cy] = size(XC2mean);
               if cx ~= cy
                    fprintf('m_MahalanobisD :: Error(0) :: Mahalanobis distanse calculation :: It is required the input to have the same number of columns\n');
                    error('m_MahalanobisD :: Error(1)');
                    return;
               end
               if rx < cx
                    fprintf('m_MahalanobisD :: Error(2) :: Mahalanobis distanse calculation ::The number of rows of X must exceed the number of columns\n');
                    error('m_MahalanobisD :: Error(2)');
                    return;
               end
               m = mean(XC1mean);
               M = m(ones(ry,1),:);
               C = XC1mean - m(ones(rx,1),:);
               [Q,R] = qr(C,0);
               if ( abs(det(R)) < err_ZeroDet )
                   if m_Display_PARAMS.show_flag
                        fprintf('m_MahalanobisD :: Warning(1) :: For Class No. %d and Class No. %d data is not homogenous, database must be splited\n',m_classIndex,m_classIndex2);
                   end;
                   d = 0; %%% data is not homogenous, not came from the same classes
                          %%% we made a supposition that the problem is complex and task must be devide
               else
                   ri = R'\(XC2mean-M)';
                   d = sum(ri.*ri)'*(rx-1); %%% Mahalanobis distance
               end;
               if ( d < 0 )
                    error('m_MahalanobisD :: Error(3) :: Mahalanobis distance is negative. Normalization of complexity ration can done correctly!');
                    return;
               else
                    m_indexFDR = m_indexFDR+1;
                    m_multiClassFDRVector(m_indexFDR) = d;
               end;
        end; %%% end of if_construction              
    end;  %%% m_classIndex2
end;  %%% m_classIndex
%%%% Normalazing with
m_ratio = min(m_multiClassFDRVector)/(1+min(m_multiClassFDRVector)); %%% Choosing the most complex indicator over multiclassComplexity*components of X(_classIndex,...)