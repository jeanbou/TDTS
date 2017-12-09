function PC = m_distrClasses(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds distribution of the classes %
% PC - array of classes' distribution              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( size(T,2) == 0 )
    error('m_distrClasses.m :: Error(1) :: Incomming array of classes is empty');
end;

if ( min(T) == max(T) ) %%%%% Case when there is only one class
    PC(1) = 1;
    return;
end;

clear PC;
k = 1;
m_sumSizeIndex = 0;
for i = min(T):max(T)
  indexes = find(T == i);
  m_sizeIndx = size(indexes,2); 
  %%% we skip this step when there is a space between classes
  if ( m_sizeIndx ~= 0 )
      PC(k) = m_sizeIndx;
      k = k + 1;
      m_sumSizeIndex = m_sumSizeIndex + m_sizeIndx; 
  end; 
end;
PC = PC/m_sumSizeIndex; 