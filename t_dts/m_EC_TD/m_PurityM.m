function m_ratio = m_PurityM(X,T,B);
% --------------------------------------------------------------
% Purity Measure, proposed by Singh based on PRISM ref: 2003Y
% Input:
% - T - classes labels, because of such data sructure, X information do bot used
% - X - database
% Output:
% - m_ration - final Purity Measure
% Copyright by Ivan Budnyk 2008 / the 1st of March
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONCTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
NORMALIZATION_PARAMETER = 0.702; % based on Singh's result
MAXIMAL_B = 31; % Maximal B in the range [0;31], propossed and verified by Singh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (B > MAXIMAL_B)
    fprintf('m_Collective_Entropy :: Note(1.1) :: B - resolution is out of tested by Sing range, B set to MAXIMAL_B\n');
    B = MAXIMAL_B+1;
end; % check of incomming paramter B

if (B < 2)
    fprintf('m_Collective_Entropy :: Note(1.2) :: B - resolution is out of tested by Sing range, B set to 2\n');
    B = 2;
end; % check of the incoming parameter B

for i_attribute = 1:size(X,1)
    XmaxAtt(i_attribute) = max(X(i_attribute,:));
    XminAtt(i_attribute) = min(X(i_attribute,:));
end;

for i_attribute = 1:size(X,1)
    array_ofBDeltas(i_attribute) = XmaxAtt(i_attribute)-XminAtt(i_attribute);  %%%% we take the minimal radius over all components (directions) of X
end;

if ( max(array_ofBDeltas) == 0 )
    fprintf('m_Collective_Entropy :: Warning(1) :: components in all direction are concentrated in one dot, problem predefined as 1(non)-complex\n');
    m_ratio = 1;
	return;
end;

clear('flag_nonZeroRadius','XmaxAtt','XminAtt');
for i_attribute = 1:size(X,1) %%% random choose of the initial centre of the 1st cluster // random choose of the X0 for all components in all directions
    array_ofBDeltas(i_attribute) = array_ofBDeltas(i_attribute)/B;    
end;

%%% BUILDING GRID
SHlKl = zeros(size(X,1),B,max(T)-min(T)+1);
for i = 1:size(X,1)
    if (array_ofBDeltas(i) ~= 0 )
        m_tmp_Bi=min(X(i,:));
        for res = 1:B
            for j = 1:size(T,2)
                if ( res == B ) %%% special sub case when we reached last cell of the attribute, in this case we need to include last border line
                    if ( (X(i,j) >= m_tmp_Bi) & (X(i,j) <= m_tmp_Bi+array_ofBDeltas(i)) )
                        SHlKl(i,res,T(j)) = SHlKl(i,res,T(j)) + 1;                 
                    end;
                else
                    if ( (X(i,j) >= m_tmp_Bi) & (X(i,j) < m_tmp_Bi+array_ofBDeltas(i)) )
                        SHlKl(i,res,T(j)) = SHlKl(i,res,T(j)) + 1;                 
                    end;
                end;
            end; %%% end by j value of
            m_tmp_Bi = m_tmp_Bi + array_ofBDeltas(i);
        end; %%% end for B-resolution
    end;
end; %%% end for component
Pil = SHlKl/size(T,2);
clear('SHlKl');
%%%% Calculating SHL
for i = 1:size(X,1)
    if (array_ofBDeltas(i) ~= 0 )
        for res = 1:B
            SHL(i,res) = 0;
            for j = min(T):max(T)
                SHL(i,res) = SHL(i,res) +  ( ( Pil(i,res,j)- 1/(max(T)-min(T)+1) )^2 );                
            end;
            SHL(i,res) = (max(T)-min(T)+1)*SHL(i,res)/(max(T)-min(T));
            SHL(i,res) = sqrt(SHL(i,res));
        end;
    end;
end;

%%%% Calculating SH
for i = 1:size(X,1)
    SH(i) = 0;
    if (array_ofBDeltas(i) ~= 0 )
        for res = 1:B
            SumN = 0;
            for j = min(T):max(T)
                SumN = SumN + Pil(i,res,j); 
            end;
            SH(i) = SH(i) + SHL(i,res) * SumN;
        end;
    end;
end;
clear('SHL','Pil');
Ec = SH/(2^(B-1));  %%% "this is to ge a larger weight for higher resolution" // the citation of article of Singh for 2004
clear('SH');
Ec = Ec/NORMALIZATION_PARAMETER;
m_ratio = min(Ec); %%% pick-up the hardest case of attribute, where the entropy maximal // means Ec - minimal, ref article of Singht