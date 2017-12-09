function ch = m_xSymb(nn);
% ----------------------------------------------------------------------------------------------
%  (c) 2005 Mariusz RYBNIK / modified by Ivan BUDNYK (c) 2008
%  Function xSymb generates combinations of plot styles required for different type of clusters
%  nn - index of incominf cluster
%  ch _ generated color and style of cluster
%  MODIFICATION: To select them in random way
% ----------------------------------------------------------------------------------------------

M_NUM_COLORS =  6;
M_NUM_STYLES = 13;


dotstyle = round(rand*(M_NUM_STYLES-1));
color = round(rand*(M_NUM_COLORS-1));

% OLD CODE / by RYBNIK MARIUSZ
%nn = rem(nn,M_NUM_COLORS*M_NUM_STYLES);
%dotstyle = fix(nn/M_NUM_COLORS);
%color = rem(nn,M_NUM_COLORS);

%%%% Selecting color
switch color
	case 0, ch = 'g';
    case 1, ch = 'b';
	case 2, ch = 'c';
	case 3, ch = 'm';
	case 4, ch = 'y';
	case 5, ch = 'k';	
end;
%%%% ADDING SELECTING DOT-STYLE
switch dotstyle
    case 0, ch = [ch  '+'];  
	case 1, ch = [ch  '.' ];
	case 2, ch = [ch  'v' ];
    case 3, ch = [ch  'x' ];
	case 4, ch = [ch  '*' ];
	case 5, ch = [ch  's' ];
	case 6, ch = [ch  'd' ]; 
    case 7, ch = [ch  'o' ]; 
	case 8, ch = [ch  '^' ];
	case 9, ch = [ch  '<' ];
	case 10, ch = [ch  '>' ];
	case 11, ch = [ch  'p' ];
	case 12, ch = [ch  'h' ];
end;