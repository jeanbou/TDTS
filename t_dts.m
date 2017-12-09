% ------------------------------------------- Main Programme -------------------------------------------%
% Tree-like Divide To Simplify                                                                          %
% Concept propossed by Prof. Kurosh Madani, MDC Chebira Nasser, Dr. Mariusz Rybnik                      %
% Current Version 2.50                                                                                  %
% Copyright(c) Ivan Budnyk 01-Jan-2009                                                                  %   
% Software Component Structure and Protoversion 2.00 Beta done by Dr. Mariusz Rybnik (c) 2005           %
%--------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT LAUNCHING OF GUI WINDOW %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this part of code we set/save/load parameters from GUI and pass to T-DTS        %
% During T-DTS development parameters has been created as constants in *.m - script  %
% So it's useful to look-up for original decription for m_ini_tdts.m.txt file        %
% Then it was reconverted to mat file. The standart of T-DTS parameters I've saved   %
% BECAUSE: 'save' command doesn't allow to append obj to file, so in interface code  %
% file M_STR_PARAMSFILENAME - each time loadding, a variable updating and then saving%
% at the end this file will consist parameters and GUI variables                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = t_dts(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global M_STR_PARAMSFILENAME
M_STR_PARAMSFILENAME = 't_dts_params.mat'   ; %%% by default name of parameters of T-DTS parameters
cd C:\MATLAB6p1\work\tdts\                    %%% by default path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF THE CONSTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0  % LAUNCH GUI
	fig = openfig(mfilename,'reuse');
	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
	if nargout > 0
		varargout{1} = fig;
	end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end;
end;

% ------------------------------ LAUNCHING T-DTS -----------------------------------------
function varargout = figure1_CreateFcn(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME
t_dts_init_script;      %%% This script consist realization of the function of GUI creating. In it's implementation this script initialize T-DTS within by default parameter

%---------------- IN THIS PART WE UPDATE T_DTS PARAMS WITH SELECTED / SWITCHED / INPUT VALUES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DECOMPOSITION PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = m_du_menu_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

%%% Pick-up DU main parameter of decomposition from interface
m_struct_du_menu = get(handles.m_du_menu,'String');
m_num_in_menu = get(handles.m_du_menu,'Value');
load(M_STR_PARAMSFILENAME);
m_DU.method = m_struct_du_menu{m_num_in_menu};
clear('m_struct_du_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

function varargout = m_du_param_menu_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

%%% Pick-up DU parameter of type of the prototype building
m_struct_du_menu = get(handles.m_du_param_menu,'String');
m_num_in_menu = get(handles.m_du_param_menu,'Value');
load(M_STR_PARAMSFILENAME);
m_DU.Dec_Param = m_struct_du_menu{m_num_in_menu};
clear('m_struct_du_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = m_pu_menu_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

%%% Pick-up PU parameter of type of the prototype building
m_struct_pu_menu = get(handles.m_pu_menu,'String');
m_num_in_menu = get(handles.m_pu_menu,'Value');
load(M_STR_PARAMSFILENAME);
m_PU.method = m_struct_pu_menu{m_num_in_menu};
clear('m_struct_pu_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

%-------------------------------- DISPLAY PROCESSING ---------------------------
function varargout = m_dsp_prcss_ch_bx_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_Value = get(handles.m_dsp_prcss_ch_bx,'Value');
max_Value = get(handles.m_dsp_prcss_ch_bx,'Max');
load(M_STR_PARAMSFILENAME);
if ( m_Value == max_Value)
    m_Display_PARAMS.show_flag = 1;       
else
    m_Display_PARAMS.show_flag = 0;                            % then checkbox is checked-take approriate action
end
%%%% MORE GENERAL SOLUTION, ALL DISPLA FLAGS ON / OR OFF
if m_Display_PARAMS.show_flag 
    m_Display_PARAMS.show_learningSlaves_batchSize = 25;       %%%% BY DEFAULT PARAMETERS CAN BE CHANGED
else                          
    m_Display_PARAMS.show_learningSlaves_batchSize = 2000;     %%%% PARAMETERS CAN BE CHANGED / SET AS SET
end;
clear('m_Value','maxValue');
%%%% Update within last check box selection
save(M_STR_PARAMSFILENAME);

% ------------------------------- SELECTING A MODE -------------------------------------
function varargout = m_tdts_run_md_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_struct_pop_menu = get(handles.m_tdts_run_md,'String');
m_num_in_menu = get(handles.m_tdts_run_md,'Value');
load(M_STR_PARAMSFILENAME);
m_str_TDTS_mode = m_struct_pop_menu{m_num_in_menu};
clear('m_struct_pop_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

% --------------------------------- INPUT NUMBER OF ITERATIONS -----------------------------------
function varargout = m_it_num_txt_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_str_it_num = get(handles.m_it_num_txt,'String');
load(M_STR_PARAMSFILENAME);
if ( isnumeric(str2num(m_str_it_num)) & (str2num(m_str_it_num) >= 1 ) )  
    m_iterationNumber = round(str2num(m_str_it_num));
else
    fprintf('T_DTS.m :: Error(0) :: Number of iterations parameter is out of range, switched to by default value 1\n');
    m_iterationNumber = 1;
end;
clear('m_str_it_num');
save(M_STR_PARAMSFILENAME);

% ---------------------------------- COMPLEXITY THRESHOLD ----------------------------------
function varargout = m_cmplx_thr_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_str_it_num = get(handles.m_cmplx_thr,'String');
load(M_STR_PARAMSFILENAME);
if ( isnumeric(str2num(m_str_it_num)) & (str2num(m_str_it_num) >= 0 ) & (str2num(m_str_it_num) <= 1 ) )      
    m_thresholdInterval = [str2num(m_str_it_num)];
    if ( str2num(m_str_it_num) == 1 )
        fprintf('T_DTS.m :: Warning(1) :: Complexity threshold is equal 1, you have requested full decomposition of DB, please be ware that it can produce for Tree following mode a Matrix Index Exceed Error for Big Database. It appears during the process of decomposition for Generalization Database\n');
    end;
else
    fprintf('T_DTS.m :: Error(1) :: Set parameter of complexity is out of range [O;1], switched to by default value 0.5\n');
    m_thresholdInterval = [.5];
end;
clear('m_str_it_num');
save(M_STR_PARAMSFILENAME);

% ----------------------------------- COMPLEXITY TYPE SELECTION ---------------------------------
function varargout = m_cmplx_lst_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_struct_pop_menu = get(handles.m_cmplx_lst,'String');
m_num_in_menu = get(handles.m_cmplx_lst,'Value');
load(M_STR_PARAMSFILENAME);
m_complexityESTList = {m_struct_pop_menu{m_num_in_menu}};
clear('m_struct_pop_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

% -----------------------------------  EDIT FILENAME OF THE RESULTS ---------------------------------
function varargout = m_rslt_file_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_tmp_structure = get(handles.m_rslt_file,'String');
load(M_STR_PARAMSFILENAME);
m_str_outputFileName = m_tmp_structure;
clear('m_tmp_structure');
save(M_STR_PARAMSFILENAME);

% ------------------------------------- EDIT INPUT DB FILENAME -------------------------------
function varargout = m_db_filename_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_tmp_structure = get(handles.m_db_filename,'String');
load(M_STR_PARAMSFILENAME);
m_str_dataFileName = m_tmp_structure;
clear('m_tmp_structure');
save(M_STR_PARAMSFILENAME);

% ------------------------------------ FLAG OF RANDOMIZATION OF INCOMING DATA --------------------------------
function varargout = m_ran_data_ch_bx_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_Value = get(handles.m_ran_data_ch_bx,'Value');
max_Value = get(handles.m_ran_data_ch_bx,'Max');
load(M_STR_PARAMSFILENAME);
if ( m_Value == max_Value)
    m_randomizeIncomingData_flag = 1;
else
    m_randomizeIncomingData_flag = 0;
end;
clear('m_Value','maxValue');
%%%% Update within last check box selection
save(M_STR_PARAMSFILENAME);

%---------------------------------------- SLECTING A MODE OF LEARNING DB EXTRACTION ----------------------------
function varargout = m_lrn_db_extr_mode_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_struct_pop_menu = get(handles.m_lrn_db_extr_mode,'String');
m_num_in_menu = get(handles.m_lrn_db_extr_mode,'Value');
load(M_STR_PARAMSFILENAME);
m_LearningDBCreating.mode = m_struct_pop_menu{m_num_in_menu};
clear('m_struct_pop_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

% ----------------------------------------- OBTAINING SIZE in % OF LEARNING DATABASE --------------------------
function varargout = m_lrn_db_rate_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_str_it_num = get(handles.m_lrn_db_rate,'String');
load(M_STR_PARAMSFILENAME);
if ( isnumeric(str2num(m_str_it_num)) & (str2num(m_str_it_num) > 0 ) & (str2num(m_str_it_num) <= 100 ) )  
    m_LearningDBCreating.percentage = str2num(m_str_it_num);
else
    fprintf('T_DTS.m :: Error(18) :: Set percentage of the size of learning database is out of range (0;100], default value automatically selected\n');
    m_LearningDBCreating.percentage = 50;
end;
clear('m_str_it_num');
save(M_STR_PARAMSFILENAME);

% ------------------------- SELECTING DATA NORMALAZING/CONVERTING  -------------------------------------------
function varargout = m_data_conv_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_struct_pop_menu = get(handles.m_data_conv,'String');
m_num_in_menu = get(handles.m_data_conv,'Value');
load(M_STR_PARAMSFILENAME);
m_str_Type_ofNormalizing = m_struct_pop_menu{m_num_in_menu};
clear('m_struct_pop_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

% --------------------------- PCA/T CHECK BOX -----------------------------------------
function varargout = m_pca_ch_bx_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_Value = get(handles.m_pca_ch_bx,'Value');
max_Value = get(handles.m_pca_ch_bx,'Max');
load(M_STR_PARAMSFILENAME);
if ( m_Value == max_Value)
    m_PCA.flag = 1;
else
    m_PCA.flag = 0;
end;
clear('m_Value','maxValue');
%%%% Update within last check box selection
save(M_STR_PARAMSFILENAME);

% ---------------------------------- INPUT THE RATION OF COMPONENT ELEMINATION ----------------------------------------
function varargout = m_pca_prc_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_str_it_num = get(handles.m_pca_prc,'String');
load(M_STR_PARAMSFILENAME);
if ( isnumeric(str2num(m_str_it_num)) & (str2num(m_str_it_num) > 0 ) & (str2num(m_str_it_num) < 100 ) )  
    m_PCA.per_elmnt_cmps = str2num(m_str_it_num);
else
    fprintf('T_DTS.m :: Error(3) :: Set percentage of the component elemination of PCT is out of range (0;1), default value 1.0 automatically set\n');
    m_PCA.per_elmnt_cmps = 1.0;
end;
clear('m_str_it_num');
save(M_STR_PARAMSFILENAME);

% -------------------------- SELECTING DISTANCE TYPE ------------------------------------------
function varargout = m_dist_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_struct_pop_menu = get(handles.m_dist,'String');
m_num_in_menu = get(handles.m_dist,'Value');
load(M_STR_PARAMSFILENAME);
m_str_Type_ofDistance = m_struct_pop_menu{m_num_in_menu};
clear('m_struct_pop_menu','m_num_in_menu');
save(M_STR_PARAMSFILENAME);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION OF MENU IMPLEMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function varargout = m_config_Callback(h, eventdata, handles, varargin)
%%%%%%% NOTHING TO DO HERE JUST FOR SELECTION OF Sub-menu of Configuration

% --------------------------- CALL DECOMPOSITION METHODS FOR MODIFICATION -----------------------------------------
function varargout = m_du_conf_Callback(h, eventdata, handles, varargin)
edit m_getWPrototype;

% ----------------------------- CALL FOR EDIT DU Parameters ---------------------------------------
function varargout = m_du_params_conf_Callback(h, eventdata, handles, varargin)
edit m_getDecMethod_params;

% --------------------------- CALL PROCESSIN METHODS FOR MODIFICATION -------------------------------------------
function varargout = m_pu_conf_Callback(h, eventdata, handles, varargin)
edit m_build_sANN;
edit m_run_sANN;

% -------------------------------CALL FOR EDIT PARAMS OF CONF METOD -------------------------------------
function varargout = m_pu_params_conf_Callback(h, eventdata, handles, varargin)
edit m_getProcMethod_params;

% ---------------------------- CONFIGURE MULTI RUN SCRIPT FOR COMPLEX TESTING ----------------------------------------
function varargout = m_multi_run_conf_Callback(h, eventdata, handles, varargin)
edit m_multi_run_config_scrpt;

% ------------------------------ CONFIG OPTIMIZATION FUNCTION --------------------------------------
function varargout = m_opt_fun_conf_Callback(h, eventdata, handles, varargin)
edit m_optimizationFunction;

% ------------------------------- CONCATENATE RESULTS -------------------------------------
function varargout = m_concat_output_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

fprintf('T_DTS.m :: Warning(2) :: Please verify settings before processing:\n- Files are numerated properly: beginning from 1, without spaces including maximal\n- Output data file name field consists only base-name without numeration\nPress any key to continue ...\n');
pause;
load(M_STR_PARAMSFILENAME);
m_flag_existingOutputFile = 1;  %%% temporary flag required for the searching of maximal file index to concatenate
m_CONST_NUMFILES = 1;
while m_flag_existingOutputFile
    if  ( exist( strcat( m_str_outputFileName,int2str(m_CONST_NUMFILES),'.mat' ) ) == 2 ) 
        m_CONST_NUMFILES = m_CONST_NUMFILES + 1;
    else        
        m_flag_existingOutputFile = 0;  
    end;
end;
clear('m_flag_existingOutputFile','m_res','m_res_stat');

load( strcat( m_str_outputFileName,int2str(1),'.mat' ) ); %%% uploadfirst file
m_res = m_results;              
m_res_stat = m_results_stat;    
clear(M_CONST_STR_ALLRESULTSFILENAME,M_CONST_STR_STATRESULTSFILENAME);
delete(strcat( m_str_outputFileName,int2str(1),'.mat' ));

for f = 2 : m_CONST_NUMFILES-1
    load( strcat( m_str_outputFileName,int2str(f),'.mat' ) ); %%% upload next file
    m_res_SZ_i = size(m_results,1);
    m_res_SZ_j = size(m_results,2);
    for j = 1 : size(m_results,2)
       if ~isempty(m_results_stat{j})
            for k = 1 : size(m_results,3)
                m_res(m_res_SZ_i,j,k) = m_results(1,j,k);
            end;
            m_res_stat{m_res_SZ_i,j} = m_results_stat{1,j};
       end;
    end;
    clear(M_CONST_STR_ALLRESULTSFILENAME,M_CONST_STR_STATRESULTSFILENAME);
    delete(strcat( m_str_outputFileName,int2str(f),'.mat' ));
end;
m_results_stat = m_res_stat;
m_results = m_res;
clear('m_res_stat','m_res','j','k','f');
save(strcat(m_str_outputFileName,'.mat'),M_CONST_STR_ALLRESULTSFILENAME,M_CONST_STR_STATRESULTSFILENAME,'POUT','WOUT','TREE','MAXTREELEVEL');
fprintf('T_DTS.m :: Concatenation has been finished, please be informed that the other data-results belong to last iteration of T-DTS\n');

% --------------------------------------------------------------------
function varargout = m_set_const_Callback(h, eventdata, handles, varargin)
%%%%%%% NOTHING TO DO HERE JUST FOR SELECTION OF Sub-menu of Conctants configuration
% --------------------------- SET CONSTANT OF SIZE ---------------------------
function varargout = m_undec_db_sz_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter accuracy of the threshold:'};
m_title   = 'Set Accuracy of the Threshold';
m_lines= 1;
load(M_STR_PARAMSFILENAME);
m_tmp_method = m_getDecMethod_params(m_DU.method);
m_def     = {num2str(min_Size)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) & (str2num(m_answer{1}) < 1 ) )
    min_Size = str2num(m_answer{1});
else
    fprintf('T_DTS.m :: Error(4) :: Set Accuracy of the Threshold is out of range (0;1). It is set equal to 0.01\n');
    min_Size = 0.01;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer','m_tmp_method');
save(M_STR_PARAMSFILENAME);

% ------------------------- SET ZERO Epsilon neighborhood -------------------------------------------
function varargout = m_eps_nhd_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter epsilon neighborhood of sero value:'};
m_title   = 'Set Epsilon Neighborhood';
m_lines= 1;
%%%% This call of procedure need to set zero meaning value for comparicing real value close to  zero
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(zeroMeaningValue)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) & (str2num(m_answer{1}) < 1 ) )
    zeroMeaningValue = str2num(m_answer{1});
else
    fprintf('T_DTS.m :: Error(5) :: Set Constant Epsilon Neighborhood of Zero is out of (0;1). This parameter set equal to 1e-012\n');
    zeroMeaningValue = 1e-012;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% ------------------------------- SET SPLITTING RATIO FOR SELFTUNING MODE-------------------------------------
function varargout = m_splt_ration_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter splitting coefficient value:'};
m_title   = 'Set Splitting Coefficient';
m_lines= 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(M_SP_COEF)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) & (str2num(m_answer{1}) < 1 ) )
    M_SP_COEF = str2num(m_answer{1});
else
    fprintf('T_DTS.m :: Error(14) :: M_SP_COEF is out of (0;1). This parameter set equal to golden ratio\n');
    M_SP_COEF = 0.61803398875; %%% golden ratio
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% -------------------------------- SET MAXIMAL THRESHOLD ------------------------------------
function varargout = m_set_max_tr_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt = {'Enter Alfa Coefficient for PDF:'};
m_title  = 'Set Alfa Coefficient';
m_lines = 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(M_CONS_MIN_THRESHOLD)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) & (str2num(m_answer{1}) < 1 ) )
    M_CONS_MIN_THRESHOLD = str2num(m_answer{1});
else
    fprintf('T_DTS.m :: Error(15) :: Alfa Coefficient is out of (0;1). This parameter set equal to 0.1\n');
    M_CONS_MIN_THRESHOLD = 0.1; %%% by default
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% --------------------------------------------------------------------
function varargout = m_es_options_Callback(h, eventdata, handles, varargin)
%%%%%%% NOTHING TO DO HERE JUST FOR SELECTION OF SIB-menu of Estimating complexity params
% --------------------------------------------------------------------
function varargout = m_ec_sse_goal_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter Distance Type'};
m_title   = 'Type of Distance for ZISC Complexity Estimator';
m_lines= 1;
load(M_STR_PARAMSFILENAME);
m_def     = {m_RBFParams4ECTD.dstn_type};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
m_RBFParams4ECTD.dstn_type = m_answer{1};
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% --------------------------------------------------------------------
function varargout = m_ec_max_data_size_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter MIF for ZISC Complexity Estimator:'};
m_title   = 'MIF';
m_lines   = 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(m_RBFParams4ECTD.MIF)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) >= 0 ) )
    m_RBFParams4ECTD.MIF = round(str2num(m_answer{1}));
else
    fprintf('T_DTS.m :: Error(7) :: Set MIF for ZISC Complexity Estimator is out of range\n');
    m_RBFParams4ECTD.MIF = 0;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% --------------------------------------------------------------------
function varargout = m_ec_it_num_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter Iteration Number for ZISC Complexity Estimator:'};
m_title   = 'Iteration Number';
m_lines   = 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(m_RBFParams4ECTD.n_iteration)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) )
    m_RBFParams4ECTD.n_iteration = round(str2num(m_answer{1}));
else
    fprintf('T_DTS.m :: Error(8) :: Set Iteration Number for ZISC Complexity Estimator is out of range\n');
    m_RBFParams4ECTD.n_iteration = 1;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% --------------------------------------------------------------------
function varargout = m_ec_pol_pow_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter Polynome Power for ZISC Complexity Estimator:'};
m_title   = 'Polynome Power';
m_lines   = 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(m_RBFParams4ECTD.m_pol_pow)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) > 0 ) )
    m_RBFParams4ECTD.m_pol_pow = round(str2num(m_answer{1}));
else
    fprintf('T_DTS.m :: Error(9) :: Set SPolynome Power for ZISC Complexity Estimator is out of range\n');
    m_RBFParams4ECTD.m_pol_pow = 4;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% --------------------------------------------------------------------
function varargout = m_ec_res_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

m_prompt  = {'Enter Resolution Parameter for Complexity Estimating Procedure:'};
m_title   = 'Set Resolution';
m_lines= 1;
load(M_STR_PARAMSFILENAME);
m_def     = {num2str(resolution_B)};
m_answer  = inputdlg(m_prompt,m_title,m_lines,m_def);
if ( isnumeric(str2num(m_answer{1})) & (str2num(m_answer{1}) >= 2 ) )
    resolution_B = round(str2num(m_answer{1}));
else
    fprintf('T_DTS.m :: Error(16) :: Set Resolution is less than two. This parameter set equal to 2\n');
    resolution_B = 2;
end;
clear('m_prompt','m_title','m_lines','m_def','m_answer');
save(M_STR_PARAMSFILENAME);

% ---------------------- MENU ANALYSIS ---------------
function varargout = m_analysis_menu_Callback(h, eventdata, handles, varargin)
%%%%%%% NOTHING TO DO HERE JUST FOR SELECTION OF Sub-menu Analysis

% ---------------------- BUILD PDF of COMPLEXITY OF MAXIMAL DECOMPOSITION OVER THRESHOLD --------------
function varargout = m_pdf_complexity_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
load(M_CONST_OPTIMUM_STRFILENAME,'m_histXInterval','m_pdfCML');

axes(handles.m_tdts_axes);
set(handles.m_tdts_axes_title,'String','Chart :: PDF of Complexity');
plot(m_histXInterval,m_pdfCML,':*b');
axis tight;
grid on;
clear('m_histXInterval','m_pdfCML');

% ------------------ CALCULATE ENTROPY OF PDF OF MAXIMAL DECOMPOSITION -------------------------------
function varargout = m_post_reg_lrn_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
load(M_CONST_OPTIMUM_STRFILENAME,'m_pdfCML');

for i_temp = 1 : size(m_pdfCML,2)
    if ( m_pdfCML(i_temp) ~= 0 )
        m_entropy(i_temp) = - m_pdfCML(i_temp)*log(m_pdfCML(i_temp)); 
    else
        m_entropy(i_temp) = 0;
    end;
end;
fprintf('* T_DTS.m :: Entropy of PDF of Max-decomposition : %f\n',sum(m_entropy));
clear('m_pdfCML','m_entropy','i_temp');

% -----------------------------The same but for generalizing ---------------------------------------
function varargout = m_post_reg_gen_Callback(h, eventdata, handles, varargin)
%global M_STR_PARAMSFILENAME

fprintf('WARNING(3) : POST REGRESSION ANALYZIS IS NOT REALIZED IN t_dts.m! IS IT REQUIRED?\n');

% --------------------------- PRINT OUT THE OPTIMAL SOLUTION -----------------------------------------
function varargout = m_prnt_fnd_opt_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
load(m_str_outputFileName,'m_results_stat');

% Searh for optimal in data file, preprocess info about threshold and print optimal threshold out out
if (size(m_results_stat,1) == 1 )
    m_opt = m_results_stat{1,size(m_results_stat,2)}.optimum;
    m_opt_Thresh = m_results_stat{1,size(m_results_stat,2)}.threshold;
    m_opt_avr_LR = m_results_stat{1,size(m_results_stat,2)}.avr_LR;
    m_opt_std_LR = m_results_stat{1,size(m_results_stat,2)}.std_LR;
    m_opt_avr_GR = m_results_stat{1,size(m_results_stat,2)}.avr_GR;
    m_opt_std_GR = m_results_stat{1,size(m_results_stat,2)}.std_GR;
    m_opt_avr_ET = m_results_stat{1,size(m_results_stat,2)}.avr_ET;
    m_opt_std_ET = m_results_stat{1,size(m_results_stat,2)}.std_ET;
    for i = 1 : (size(m_results_stat,2)-1)
        if ~isempty(m_results_stat{1,i})
            if (m_results_stat{1,i}.optimum < m_opt)
                m_opt = m_results_stat{1,i}.optimum;
                m_opt_Thresh = m_results_stat{1,i}.threshold;
                m_opt_avr_LR = m_results_stat{1,i}.avr_LR;
                m_opt_std_LR = m_results_stat{1,i}.std_LR;
                m_opt_avr_GR = m_results_stat{1,i}.avr_GR;
                m_opt_std_GR = m_results_stat{1,i}.std_GR;
                m_opt_avr_ET = m_results_stat{1,i}.avr_ET;
                m_opt_std_ET = m_results_stat{1,i}.std_ET;                
            end;
        end;        
    end;
    fprintf('T_DTS.m :: OPTIMUM\n-------------------------------------------------------------\n');
    fprintf('* Optimal TDTS executing time : %4.4f +/- %4.4f\n* Optimal Learning rate       : %4.4f +/- %4.4f\n* Optimal Generalization rate : %4.4f +/- %4.4f\n',m_opt_avr_ET,m_opt_std_ET,m_opt_avr_LR,m_opt_std_LR,m_opt_avr_GR,m_opt_std_GR);
    fprintf('-------------------------------------------------------------\n* Optimal threshold : %1.4f\n',m_opt_Thresh);
            
else
    fprintf('T_DTS.m :: Error(6) :: This code realized for only one complexity threshold, please modify parameters of multi run, and provide retesting\n');
end;

% ------------------------- DRAW OPTIMIZATION FUNCTION -------------------------------------------
function varargout = m_optm_fun_chrt_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
load(m_str_outputFileName,'m_results_stat');

%%% drawing chart of optimization function
tmp_jIndx = 1;
if (size(m_results_stat,1) == 1 )
    for i = 1 : size(m_results_stat,2)
        if ~isempty(m_results_stat{1,i})            
                m_optm(tmp_jIndx) = m_results_stat{1,i}.optimum;
                m_opt_Threshl(tmp_jIndx) = m_results_stat{1,i}.threshold;
                tmp_jIndx = tmp_jIndx + 1;
        end;        
    end;
    m_handle_opt_fig = figure('Name','Threshold analytical chart :: Function - Optimum(Threshold)','numbertitle','off');        
    plot(m_opt_Threshl,m_optm,'dr');
    %axis([0 0 0.9 38]);
    grid on;
else
    fprintf('T_DTS.m :: Error(17) :: This function realized for only one complexity threshold, please modify parameters of multi run, and provide retesting\n');
end;



% ------------------------- HELP MENU -------------------------------------------
function varargout = m_help_Callback(h, eventdata, handles, varargin)

% ----------------------------- ABOUT ---------------------------------------
function varargout = m_about_Callback(h, eventdata, handles, varargin)
helpdlg('(c) Ivan BUDNYK 2009','About T-DTS v2.50');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPLEMENTATION OF THE "RUN/GO/DO" BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALL OF TDTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = m_tdts_run_button_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

fprintf('\nT_DTS.m :: TDTS processing\n');
load(M_STR_PARAMSFILENAME);

%------------------------ End of Constant lists ---------------------
m_TDTS(m_CONST_MAX_PRC, M_CONST_CLASS_LABEL_UNDETERM,M_CONST_THRESHOLD_DELTA,M_CONST_OPTIMUM_STRFILENAME,M_SP_COEF,M_CONS_MIN_THRESHOLD,m_str_dataFileName,m_Display_PARAMS,m_randomizeIncomingData_flag,min_Size,zeroMeaningValue,m_RBFParams4ECTD,m_str_Type_ofDistance,m_str_Type_ofNormalizing,m_PCA,m_DU,m_PU,m_str_TDTS_mode,m_str_outputFileName,m_iterationNumber,m_thresholdInterval,m_complexityESTList,resolution_B,m_LearningDBCreating);
fprintf('T_DTS.m :: TDTS has finished processing\nPlease remove/delete %s.MAT file or\nPlease input new name of an output file for the results of the testing',m_str_outputFileName);

% ------------------------- THIS BUTTON PRINTS OUT COMPLEXITY of the solid, not decomposed DB
function varargout = m_db_complexity_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
fprintf('\nm_print_DB_Coplexity :: Complexity ::\n-------------------------------------------------------------\n');
print_DB_Complexity(m_str_dataFileName,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS,m_randomizeIncomingData_flag, m_LearningDBCreating, m_CONST_MAX_PRC, m_str_Type_ofNormalizing, m_PCA);
fprintf('------------------------------------------------------------\n');

% ------------------------- PRINTS SUMMARY RESULT -------------------------------------------
function varargout = m_rslt_btn_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
print_Results_scrpt;

% --------------------------BUILD A CHART FUCNTIONS of RESULTS FROM THRESHOLD ------------------------------------------
function varargout = m_res_thresh_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
m_struct_handles.axes = handles.m_tdts_axes;
m_struct_handles.title = handles.m_tdts_axes_title;

%%%% Load results from output file
load(m_str_outputFileName,'m_results','m_results_stat');
m_new_jInd = 1;
for j_thrsh = 1 : size(m_results_stat,2)
    %%% building arrays for plotting using fist component
    if ~isempty(m_results_stat{1,j_thrsh})
            m_Tha(m_new_jInd)  = m_results_stat{1,j_thrsh}.threshold;
            avr_GR(m_new_jInd) = m_results_stat{1,j_thrsh}.avr_GR;
            std_GR(m_new_jInd) = m_results_stat{1,j_thrsh}.std_GR;
            avr_LR(m_new_jInd) = m_results_stat{1,j_thrsh}.avr_LR;
            std_LR(m_new_jInd) = m_results_stat{1,j_thrsh}.std_LR;
            m_new_jInd = m_new_jInd + 1;;
    end;
end;
%%% Drawing a chart
cla;
axes(m_struct_handles.axes);	
axis auto;
set(m_struct_handles.title,'String','Function - Results(Threshold)');
%%% RED _ Learning, BLACK _ GENERALIZATION
plot(m_Tha,avr_GR,'db',m_Tha,avr_GR+std_GR,'*b',m_Tha,avr_GR-std_GR,'*b',m_Tha,avr_LR,'dr',m_Tha,avr_LR+std_LR,'xr',m_Tha,avr_LR-std_LR,'xr');
grid on;
clear('m_Tha','avr_GR','std_GR','avr_LR','std_LR','m_results','m_results_stat','m_new_jInd','j_thrsh');

% ------------------------- PRINT COMPLEXITY OF EACH CLUSTER ------------------------------------------
function varargout = m_clust_cmplx_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
if strcmp(m_str_TDTS_mode,'Uni')
    %%%% Load CLUSTERS from result file
    load(m_str_outputFileName,'COUT','POUT');
    fprintf('\nm_print_SubDBs_Coplexity :: %s ::\n-------------------------------------------------------------\n',m_complexityESTList{1});
    %%%% Creating temporary complexity structure for m_EC_TD call function
    m_tmp_Complexity_PARAMS.Method_name = m_complexityESTList{1};
    m_tmp_Complexity_PARAMS.Threshold = 1;
    for temp_i = 1 : size (COUT,2)
        [di,cm] = m_EC_TD(POUT{temp_i},COUT{temp_i},1,m_tmp_Complexity_PARAMS,zeroMeaningValue,m_RBFParams4ECTD,resolution_B,m_Display_PARAMS);
        fprintf('Cluster No.: %d Size : %d Complexity : %1.6f\n',temp_i,size(COUT{temp_i},2),cm);
    end;
    clear('m_tmp_Complexity_PARAMS','temp_i','di','cm');
    fprintf('------------------------------------------------------------\n');    
    clear('POUT','COUT');
else
    fprintf('\nT_DTS.m :: Error(10) :: This function can be used only for Uni mode of T-DTS run when one complexity estimating function is selected\n');
end;
save(M_STR_PARAMSFILENAME);

% ------------------------- 2D Clusters ploting -------------------------------------------
function varargout = m_2d_graph_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);

M_CONST_FONT_SIZE =   9;
M_DELTA_CENT_PROT = 0.1;

if strcmp(m_str_TDTS_mode,'Uni')
    m_struct_handles.axes = handles.m_tdts_axes;
    m_struct_handles.title = handles.m_tdts_axes_title;
    %%%% Load CLUSTERS from result file
    load(m_str_outputFileName,'POUT','WOUT');
    % Plotting in 2D clusters after decomposition
    if (size(POUT{1},1) == 2 )
        cla;
        axes(m_struct_handles.axes);
        %%% fit axis up to the ranges of 2D data
        [minX1, maxX1, minX2, maxX2] = m_getRange2D(POUT);
        axis([minX1, maxX1, minX2, maxX2]);
        axis on;
        grid off;
        set(m_struct_handles.title,'String','Clustering');
        hold on;
        for i = 1 : size(POUT,2)
            if ~isempty(POUT{i})  
            %%% Plot first and second component
            plot(POUT{i}(1,:),POUT{i}(2,:),m_xSymb(i-1));
            plot(WOUT(1,i),WOUT(2,i),'ro');
            set(text(WOUT(1,i)+M_DELTA_CENT_PROT,WOUT(2,i)-M_DELTA_CENT_PROT,sprintf('%d',i)),'FontSize',M_CONST_FONT_SIZE,'Color','red');
        end;    
    end;
    hold off;
else
    fprintf('T_DTS.m :: Error(2) :: This option is developed for two component database only\n');
end;
    clear('POUT','WOUT');
else
    fprintf('\nT_DTS.m :: Error(11) :: m_2d_graph.m :: This option is applicable only for Uni mode results, it does not provide avarage chart, because there is no avarage decomposition output database\n');
end;
save(M_STR_PARAMSFILENAME);

% ---------------------------- 2D TREE BUILDING ----------------------------------------
function varargout = m_2d_tree_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
if ~strcmp(m_DU.method,'none')
    m_Position = get(handles.m_tdts_axes,'Position');
    %%%% Load CLUSTERS from result file
    load(m_str_outputFileName,'TREE','MAXTREELEVEL');
    %%%% Initializing parameters for recursive 2D tree drawing
    Xcur = m_Position(3)/2;
    Ycur = m_Position(4);
    deltaY = m_Position(4)/MAXTREELEVEL;
    cla; %%%% Preparing AXES
    axes(handles.m_tdts_axes);    
	axis([0 m_Position(3) 0 m_Position(4)]);
    set(handles.m_tdts_axes_title,'String','2D Tree');
    %%%% Recursive drawing of the treee
    hold on;
    grid off;
    axis off;
    m_2d_tree(Xcur,Ycur,Xcur,Ycur,1,TREE,size(TREE.p,1),m_Position(3),deltaY);
    hold off;
    clear('TREE','MAXTREELEVEL','m_Position','Xcur','Ycur','deltaY');
else
    fprintf('\nT_DTS.m :: Error(12) :: m_2d_tree.m :: This option is unapplicable for selected mode, please check your DU and Run mode parameters\n');
end;
save(M_STR_PARAMSFILENAME);

% ----------------------- BUILDING TREE IN 3rd dimention over 2D clusters ------------------------
function varargout = m_3d_graph_Callback(h, eventdata, handles, varargin)
global M_STR_PARAMSFILENAME

load(M_STR_PARAMSFILENAME);
if ~strcmp(m_DU.method,'none')
    m_struct_handles.axes = handles.m_tdts_axes;
    m_struct_handles.title = handles.m_tdts_axes_title;
    %%%% Load CLUSTERS from result file
    load(m_str_outputFileName,'POUT','WOUT','TREE','MAXTREELEVEL');
    m_3d_graph_clusters(m_struct_handles,POUT,WOUT,MAXTREELEVEL,TREE);    
    clear('POUT','WOUT','TREE','MAXTREELEVEL');
else
    fprintf('\nT_DTS.m :: Error(13) :: This option is unapplicable for selected mode, please check your DU and Run mode parameters\n');
end;
save(M_STR_PARAMSFILENAME);