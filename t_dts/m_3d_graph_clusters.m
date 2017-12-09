function m_3d_graph_clusters(m_axes_handles,POUT,WOUT,Maxlevel,TREE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots in 2D clusters after decomposition
% then over them builds in 3D tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSIST ONE SUB FUNCTION

M_CONST_FONT_SIZE  =   9; %%% Size of the font for numeration of the clusters
M_DELTA_CENT_PROT  = 0.1; %%% Delta of coordinate to shift for text plotting
M_CONST_DELTA_Z    = 0.4; %%% Delta of Z-level for drawing tree in 3D, this coeficient is empirically selected, because of the ancle of 3D in 2D axes
M_CONST_DELTA_ZRT  =-0.0; %%% This constanta shifts the root of the 3D tree, required because of 3D is in 2D and some low_level tree nodes laying down and it's hard to recognized when they've mixed with data-patterns

if (size(POUT{1},1) == 2 )    
    %%%% Find maximal and minimal of two components
    cla;
    axes(m_axes_handles.axes);   
    set(m_axes_handles.title,'String','Clustering in 3D');
    grid on;
    %%%% get max and min range
    [minX1 maxX1 minX2 maxX2] = m_getRange2D(POUT);
    axis([minX1 maxX1 minX2 maxX2 0 Maxlevel]);
    axis off;
    hold on;
    for i = 1 : size(POUT,2)
      if ~isempty(POUT{i})
        %%% Creating a 3rd dimention zero flow for 3D plotting
        clear m_zeroZVector;
        m_zeroZVector(1:size(POUT{i},2)) = 0;
        %%% Plot first and second component
        plot3(POUT{i}(1,:),POUT{i}(2,:),m_zeroZVector,m_xSymb(i-1));
        plot3(WOUT(1,i),WOUT(2,i),m_zeroZVector,'ro');
        set(text(WOUT(1,i)+M_DELTA_CENT_PROT, WOUT(2,i)-M_DELTA_CENT_PROT, sprintf('%d',i)),'FontSize',M_CONST_FONT_SIZE,'Color','red');
      end;    
    end;
    %%%% Initializing parameters of the root /problamatic because 3D in 2D / of the recursive Tree drawing in 3D
    Xcur = maxX1;
    Ycur = maxX2;
    Zcur = Maxlevel+M_CONST_DELTA_ZRT;
    
    %deltaY = m_Position(4)/MAXTREELEVEL;    
    %%%% Recursive drawing of the treee
    m_3d_tree(Xcur,Ycur,Zcur,Xcur,Ycur,Zcur,TREE,size(TREE.p,1),M_CONST_DELTA_Z);
    hold off;    
else
    fprintf('m_3d_graph_clusters :: Error(1) :: This option is developed for two component database only\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m_3d_tree(Xcur,Ycur,Zcur,Xprev,Yprev,Zprev,TREE,NP,M_CONST_DELZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function recursevily builds a tree of decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot3([round(Xprev) round(Xcur)],[round(Yprev) round(Ycur)],[round(Zprev) round(Zcur)],'-*r');
Xprev = Xcur;
Yprev = Ycur;
Zprev = Zcur;
for i = 1 : NP
    Xcur = TREE.p(i,1);
    Ycur = TREE.p(i,2);
    Zcur = Zcur-M_CONST_DELZ;
    if ( ~isstruct(TREE.nnf{i}) | (NP == 1) ) %%% Now it's not tree, but a node
        plot3([round(Xprev) round(Xcur)],[round(Yprev) round(Ycur)],[round(Zprev) 0],'-sb');        
    else %%% this node consists tree, so go on
        m_3d_tree(Xcur,Ycur,Zcur,Xprev,Yprev,Zprev,TREE.nnf{i},size(TREE.nnf{i}.p,1),M_CONST_DELZ);
    end;
end;