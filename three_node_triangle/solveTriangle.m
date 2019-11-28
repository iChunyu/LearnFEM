% d = solveTriangle(K,F,cons) introduces constraint and solves FEM
%     K   --- total stiffness matrix
%     F   --- force vector on nodes
%    cons --- displacement constraint, constains [dof,u]

%% main
function d = solveTriangle(K,F,cons)
    for n = 1:length(cons(:,1))
        idx = cons(n,1);
        K(idx,:) = 0;
        K(:,idx) = 0;
        K(idx,idx) = 1;
        F = F-K(:,idx)*cons(2);
        F(idx) = cons(n,2);
    end
    d = K\F;
end