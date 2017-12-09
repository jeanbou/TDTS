function m_ratio = m_bayes_error(X,T);

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

clear('m_sortX');
for m_Class1 = min(T) : max(T)
    indexCL_tmp = find(m_Class1 == T);
    indexCLk{m_Class1} = indexCL_tmp;
    pc(m_Class1) = size(indexCL_tmp,2)/size(T,2);
    pxc{m_Class1} = zeros(1,size(pureXsort,2));  %%% preparing building p(x|c1), p(x|c2) and etc distribution
end;

clear('indexCL_tmp');
for m_Class1 = min(T) : max(T) %%% for all class distribution
    for i = 1 : size(indexCLk{m_Class1},2) %%% for every ptototype that belongs to class m_Class1       
        for j = 1 : size(pureXsort,2) %%% building pxc of class k probability, finding the index
           if isequal( X(:,indexCLk{m_Class1}(i)), pureXsort(:,j) )
              pxc{m_Class1}(j) = pxc{m_Class1}(j) + 1;
              break;
           end; 
        end;        
    end; %%% end for the propototypes of the certain class    
    pxc{m_Class1} = pxc{m_Class1}/size(indexCLk{m_Class1},2); %%% no distribution is correct        
end; %%% end for the class distribution

%%% calculating PCkX
clear('pckx','i','indexCLk');
for j = 1 : size(pureXsort,2)   
    m_sum = 0;
    for m_Class1 = min(T) : max(T)
        m_sum = m_sum + pxc{m_Class1}(j)*pc(m_Class1);        
    end;    
    for m_Class1 = min(T) : max(T)
        pckx(m_Class1,j) = pxc{m_Class1}(j)*pc(m_Class1)/m_sum;
    end;    
end;

%%% calculating bayes error
clear('pxc','pc','m_sum','m_Class1');
bayes_error = 0;
for j = 1 : size(pureXsort,2)   
    bayes_error = bayes_error + (1-max(pckx(:,j)))*px(j);
end;

m_ratio = 1-bayes_error;
clear('bayes_error','j','pckx','px');
return;