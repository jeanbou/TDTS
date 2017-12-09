function m_2d_tree(Xcur,Ycur,Xprev,Yprev,level,TREE,NP,maxWd,deltaY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function recursevily builds a tree of decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Consanta required for little Y coordinate shifting
M_CONST_DELTA_DELTA = 0.2;
%%% Ploting tree: list of last node

if (NP == 1) %%% Plot the end node
    plot([round(Xprev) round(Xcur)],[round(Yprev-M_CONST_DELTA_DELTA) round(Ycur-M_CONST_DELTA_DELTA)],'-sb');
else %%% Plot the sub-tree-node and go on
    plot([round(Xprev) round(Xcur)],[round(Yprev-M_CONST_DELTA_DELTA) round(Ycur-M_CONST_DELTA_DELTA)],'-*r');
    for i = 1 : NP
        Xprev = Xcur;
        Yprev = Ycur;
        %%% Coeficient of shifting, found empirically, however there should be explanation        
        shift_coef = 2/(NP-1);
        %%% Shifting X coordinate from the root X / required for drawing
        Xleft = Xcur - maxWd/(shift_coef*(NP^level+1));
        if ~isstruct(TREE.nnf{i}) %%% Now it's not tree, but a node
            m_2d_tree(Xleft+(i-1)*maxWd/(NP^level+1),Ycur-deltaY,Xprev,Yprev,level+1,TREE.nnf{i},1,maxWd,deltaY);
        else %%% this node consists tree, so go on
            m_2d_tree(Xleft+(i-1)*maxWd/(NP^level+1),Ycur-deltaY,Xprev,Yprev,level+1,TREE.nnf{i},size(TREE.nnf{i}.p,1),maxWd,deltaY);
        end;
    end;
end;