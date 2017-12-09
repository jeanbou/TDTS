function m_optValue = m_optimizationFunction(m_CONST_MAX_PRC, QSRate,stdQSRate,QLRate,stdQLRate,Time,stdTime,Nprot,stdNprot,m_LearningDBCreating);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My function of optimization , returns ratio of optimum. 
% Need to be minimized
%%%%%%%%%%%%%%%%%%%%%%%%% List of priorities %%%%%%%%%%%%%%%%%%%%%

M_CONST_priority1  = 1;
M_CONST_priority2  = 2;
M_CONST_priority3  = 0;
M_CONST_priority4  = 0;
M_CONST_priority5  = 0;
M_CONST_priority6  = 0;

m_optValue = ( (m_CONST_MAX_PRC-QSRate)*M_CONST_priority2 + (m_CONST_MAX_PRC-QLRate)*M_CONST_priority3 + abs(2*m_CONST_MAX_PRC-QLRate*m_LearningDBCreating.percentage/m_CONST_MAX_PRC-QSRate*(m_CONST_MAX_PRC-m_LearningDBCreating.percentage)/m_CONST_MAX_PRC-stdQLRate-stdQSRate)*M_CONST_priority1 + stdQLRate*M_CONST_priority4 + stdQSRate*M_CONST_priority5 + (Time/Nprot + stdTime/(1+stdNprot))*M_CONST_priority6  ) / (M_CONST_priority1+M_CONST_priority2+M_CONST_priority3+M_CONST_priority4+M_CONST_priority5+M_CONST_priority6);

if (m_optValue == 0)
    fprinf('m_optimizationFunction.m :: Warning(1) :: OPT Quality value is 0\n');
end