function m_ratio = m_JSD(X,T);

clear('i','j','pureXsort','m_sortX','px','m_genIndex','m_Class1','m_Class2','px');
m_sortX = sortrows(X')'; %%% 

 j = 1;
 pureXsort(:,j) = m_sortX(:,j); % create pureX -array of X-multi flat over which we will build distributions of the classes
 px(j) = 1;                     % build distribution of px(j) over pureXsort - array
 i = 2;
 while ( i <= size(m_sortX,2) ) 
     if ~isequal(pureXsort(:,j), m_sortX(:,i))
         j=j + 1;
         pureXsort(:,j) = m_sortX(:,i);
         px(j) = 1;
     else
         px(j) = px(j) + 1;
     end;
     i=i+1;
end;
%%% finishing forming general distribution px
px = px/size(T,2);
   
clear('m_genIndex');
m_genIndex = 1; %%% index that calculates pairs of class compared processing

for m_Class1 = min(T) : (max(T)-1)
    for m_Class2 = (m_Class1+1) : max(T)  %%% create all possible pair comparising

        clear('prC1','prC2','PxC1','PxC2','indexesCL1','indexesCL2','m_indexOFINCL1','m_indexOFINCL2','indexINPureX');

        PxC1 = zeros(1,size(pureXsort,2));
        PxC2 = zeros(1,size(pureXsort,2));
        indexesCL1 = find(m_Class1 == T);
        indexesCL2 = find(m_Class2 == T);
        prC1 = size(indexesCL1,2)/size(T,2);
        prC2 = size(indexesCL2,2)/size(T,2);                
        %%% building PxC1
        for m_indexOFINCL1 = 1 : size( indexesCL1,2 )
            for indexINPureX = 1 : size(pureXsort,2)
                if isequal(X(:,indexesCL1(m_indexOFINCL1)), pureXsort(:,indexINPureX))
                   PxC1(indexINPureX) = PxC1(indexINPureX) + 1;
                   break;
                end;
            end;
        end;
        %%% building PxC2
        for m_indexOFINCL2 = 1 : size( indexesCL2,2 )
            for indexINPureX = 1 : size(pureXsort,2)
                if isequal(X(:,indexesCL2(m_indexOFINCL2)), pureXsort(:,indexINPureX))
                   PxC2(indexINPureX) = PxC2(indexINPureX) + 1;
                   break;
                end;
            end;
        end;
        %%% finishing forming distribution PxC1, PxC2
        clear('tmp_counterDots');
        tmp_counterDots = sum(PxC1);
        PxC1 = PxC1 / tmp_counterDots;
        tmp_counterDots = sum(PxC2);
        PxC2 = PxC2 / tmp_counterDots;
        
        %%% calculating t_divergence for the selected pair of classes
        clear('t_divergence');
        for i = 1 : size(pureXsort,2)
            if ( (PxC2(i) ~= 0) & (PxC1(i) ~= 0) )
               t_divergence{m_genIndex,i} =  0.5*PxC1(i)*log( 2*PxC1(i)/(PxC1(i)+PxC2(i)) ) + 0.5*PxC2(i)*log( 2*PxC2(i)/(PxC1(i)+PxC2(i)) ) ;
            else
               t_divergence{m_genIndex,i} = 0;
            end;
        end;
        
        m_genIndex = m_genIndex + 1; %%% calculation for next pair of two not the same classes         
        
    end;  %%% my class Index2
end;  %%% my class Index 

clear('m_t_divergence','i','j','m_int_sum');

%%% calculating bayes error for each pair of classes
for j = 1 : size(t_divergence,1) 
    m_int_sum = 0;
    for i = 1 : size(t_divergence,2)      
        m_int_sum = m_int_sum + t_divergence{j,i};
    end;
    m_t_divergence(j) = m_int_sum;
end;

m_ratio = 1- max(m_t_divergence);
return;