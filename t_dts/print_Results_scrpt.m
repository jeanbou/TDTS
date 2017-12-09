%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script prints a summary of last results saved in file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% load reguired data for printing results out
load(strcat(m_str_outputFileName,'.mat'),'m_results','m_results_stat');
for i_cmplx = 1 : size(m_results,1)
    for j_thrsh = 1 : size(m_results,2)
        if ~isempty(m_results_stat{i_cmplx,j_thrsh})
            fprintf('\nRESULTS :: Complexity Ind.No.: %d Threshold Ind.No.: %d Iterations No. : %d\n',i_cmplx,j_thrsh,size(m_results,3));        
            fprintf('-------------------------------------------------------------\n');
            fprintf('* Average TDTS executing time : %4.4f +/- %4.4f\n* Average Number of Prototypes: %4.4f +/- %4.4f\n* Average Learning rate       : %4.4f +/- %4.4f\n* Average Generalization rate : %4.4f +/- %4.4f\n',m_results_stat{i_cmplx,j_thrsh}.avr_ET,m_results_stat{i_cmplx,j_thrsh}.std_ET,m_results_stat{i_cmplx,j_thrsh}.avr_NPs,m_results_stat{i_cmplx,j_thrsh}.std_NPs,m_results_stat{i_cmplx,j_thrsh}.avr_LR,m_results_stat{i_cmplx,j_thrsh}.std_LR,m_results_stat{i_cmplx,j_thrsh}.avr_GR,m_results_stat{i_cmplx,j_thrsh}.std_GR);
            fprintf('-------------------------------------------------------------\n');
        end;
    end;
end;
clear('m_results','m_results_stat','i_cmplx','j_thrsh');