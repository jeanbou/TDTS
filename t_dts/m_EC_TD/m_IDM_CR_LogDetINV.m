function m_ratio = m_IDM_CR_LogDetINV(X,T,err_ZeroDet);
% --------------------------------------------------------------
% Interclass Distance Measure, Criterion J2 (ref: on Fukunaga) 
% Implement for multiclass case:
% 1. Compute a ration for each combination of classes by the propossed
% 2. Choose the maximum(minimum) of the combinations
% Input:
% - X(classIndex,data_componentIndex)
% - T(classIndex) = Number/label of class
% Output:
% - m_ration
% Copyright by Ivan Budnyk 2008 / the 8th of March
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_maxClassNumber = max(T);
m_minClassNumber = min(T);
clear('m_multiClassFDRVector','m_indexFDR');
m_indexFDR = 0;
m_multiClassFDRVector(1) = 0;
for m_classIndex = m_minClassNumber:m_maxClassNumber,
	for m_classIndex2 = (m_classIndex+1):m_maxClassNumber,
       indexes1 = find(T == m_classIndex);
       indexes2 = find(T == m_classIndex2);
       if (size(indexes1) ~= 0) & (size(indexes2) ~= 0)  %%% This check allows to have holes between class Indexation, for example T = [1 3 4], this procedure is working when no class #2
               % NOTE: Each component of X consist in vector /NOTE MATRIX/ of the data, if you have a need to use X-multidimensional datacube, X(_classIndex, .... ) must be convert to datatype X(_classIndex,_dataIndex)
               clear('XC1','XC2','XC1mean','XC2mean','Xnew','GC1','GC2','GB1','GW','G');
               XC1 = X(:,indexes1);
               XC2 = X(:,indexes2);
               Xnew = XC1;
               for i=1:size(XC2,2)
                   Xnew(:,size(XC1,2)+i)=XC2(:,i);
               end; %%% Concatenation two matrix of the 2 classes only
               
               for m_ComponentNumber=1:size(XC1,1)
                   XC1mean(m_ComponentNumber) = mean(XC1(m_ComponentNumber,:));
               end;
               
               for m_ComponentNumber=1:size(XC2,1)
                   XC2mean(m_ComponentNumber) = mean(XC2(m_ComponentNumber,:));
               end;
               
               XC1mean = XC1mean';
               XC2mean = XC2mean';
                             
               %--------- Probabilities of the classes
               PC1 = size(XC1,2)/size(Xnew,2);
               PC2 = size(XC2,2)/size(Xnew,2);
               %--------- Buildning intra-matrix for Class No. 1 
               GC1 =  0;
               for m_ArrayNumber=1:size(XC1,2)
                    GC1 = GC1+(XC1(:,m_ArrayNumber)-XC1mean)*(XC1(:,m_ArrayNumber)-XC1mean)';
               end;
               %--------- Buildning intra-matrix for Class No. 2
               GC2 =  0;
               for m_ArrayNumber=1:size(XC2,2)
                    GC2 = GC2+(XC2(:,m_ArrayNumber)-XC2mean)*(XC2(:,m_ArrayNumber)-XC2mean)';
               end;
               %--------- Buildning GW intra-matrix for Class No. 1 and No. 2
               GW = PC1*GC1+PC2*GC2;
               M0 = PC1*XC1mean + PC2*XC2mean;
               %--------- Buildning GB global-matrix for Class No. 1 and No. 2
               GB1 = PC1*(XC1mean-M0)*(XC1mean-M0)'+PC2*(XC2mean-M0)*(XC2mean-M0)';
               %--------- Buildning GB global-inter/intra matrix for Class No. 1 and No. 2
               G = GW+GB1;
               
        
               if (abs(det(GW)) > err_ZeroDet)
                   tmpRatio = det(inv(GW)*G); %%% this criterion change from 1 to plus internity, 1 means simple
                   if (tmpRatio >= 1)
                        m_indexFDR = m_indexFDR + 1;
                        m_multiClassFDRVector(m_indexFDR) = log(tmpRatio); %%% Normalazing by the principle 1/(1+1). log(X) where X = [0 1], will be (0;-enternity)
                   else
                       fprintf('m_IDM_CR_LogDetINV :: Error(0.1) :: Incorrect complexity ratio, argument for logarithm is %1.4f lower than 1\n',tmpRatio);
                   end;
               else
                  fprintf('m_IDM_CR_LogDetINV :: Warning(1.1) :: Global matrix cannot be inverted, trying to use pseudoinverted matrix\n');
                  %%% in this case we use pseudo inverted matrix, approached proposed in book Karray and De Silva "Soft computing and ...", page 266 equation (5.22)
                  if (abs(det(GW'*GW)) > err_ZeroDet)
                        psINV_GW = inv(GW'*GW)*GW';
                        tmpPSINVRatio = det(inv(psINV_GW)*G);
                        if (tmpPSINVRatio >= 1)
                            m_indexFDR = m_indexFDR + 1;
                            m_multiClassFDRVector(m_indexFDR) = log(tmpPSINVRatio);
                        else
                            fprintf('m_IDM_CR_LogDetINV :: Error(0.2) :: Incorrect complexity ratio for pseudoinverted matrix, argument for logarithm is %1.4f lower than 1\n',tmpRatio);
                        end;
                  else
                      fprintf('m_IDM_CR_LogDetINV :: Warning(1.2) :: Pseudoinverted matrix cannot be inverted, problem is simple, assigning complexity ratio skipped\n');
                  end;
               end;
              
       end; %%% end of if_construction              
	end;  %%% m_classIndex2
end;  %%% m_classIndex
m_ratio = min(m_multiClassFDRVector)/(1+min(m_multiClassFDRVector)); %%% Choosing the most complex indicator over multiclassComplexity*components of X(_classIndex,...)