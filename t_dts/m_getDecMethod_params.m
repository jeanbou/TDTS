function method = m_getDecMethod_params(Dec_meth_name);
% -------------------------------------------------------------------------------------------------
% Function that determines inner parameters of each decomposition method: can be static or dynamic
% Recebtly it's statick with data-maning fuzzy abilities, it's returns method-structure
% INPUT Dec_meth_name
%   'CN' - Competitive Network
%   'SOM' - Kohonen SOM
%   'LVQ1'/'LVQ2_1' - Learning Vector Quantisizing
% INPUT complex - ratio of complexity
% --------------------------------------------------------------------------------------------------

switch Dec_meth_name
    
	case {'CNN','LVQ1','LVQ2_1'}		
		method.NoN = 2; % Totoal Number of Cluster for each splitting iteration / for GRID required method this number will be splited on number of database vector dimention
    
    case 'SOM'
        method.Grid = [1 2];
        method.NoN = 1;
        for m_tmp_ind = 1 : size(method.Grid, 2)
            method.NoN = method.NoN * method.Grid(m_tmp_ind);
        end;
        clear('method.Grid');
        
	case 'none'
		% No decomposition, no returning value
        
	otherwise 
        error('m_getDecMethod_params.m :: Error :: Not defined decomposition method name');
	end;
    
method.name = Dec_meth_name; %%% at least we determine the name of the method