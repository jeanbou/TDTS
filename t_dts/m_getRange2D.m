function [minX1, maxX1, minX2, maxX2] = m_getRange2D(POUT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function found max and min coordinates for 2D data of decomposed Database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxX1 = max(POUT{1}(1,:));
minX1 = min(POUT{1}(1,:));
maxX2 = max(POUT{1}(2,:));
minX2 = min(POUT{1}(2,:));

for i = 2 : size(POUT,2)        
  if ~isempty(POUT{i})
    if  ( maxX1 < max(POUT{i}(1,:)) )
            maxX1 = max(POUT{i}(1,:));
    end;        
    if  ( minX1 > min(POUT{i}(1,:)) )
            minX1 = min(POUT{i}(1,:));
    end;
    if  ( maxX2 < max(POUT{i}(2,:)) )
          maxX2 = max(POUT{i}(2,:));
    end;        
    if  ( minX2 > min(POUT{i}(2,:)) )
          minX2 = min(POUT{i}(2,:));
    end;        
   end;
end;