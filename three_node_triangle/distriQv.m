% Fn = distriFc(node,elem,ielem,Fe) distribute distributed force to nodes
%     node  --- coordinate of nodes in observer's frame
%     elem  --- node index of triangle elements
%     Q     --- distributed force, contains [elemIdx Qx Qy]

%% main
function Fn = distriQv(node,elem,Q)
    [Nnode,~] = size(node);
    Fn = zeros(Nnode*2,1);
    
    for n = Q(:,1)'
        x1 = node(elem(n,1),1);
        y1 = node(elem(n,1),2);
        x2 = node(elem(n,2),1);
        y2 = node(elem(n,2),2);
        x3 = node(elem(n,3),1);
        y3 = node(elem(n,3),2);
        
        A = det([1 x1 y1; 1 x2 y2; 1 x3 y3])/2;
        A = abs(A);
        Fe = Q(n,[2 3])*A;
        
        nFx = elem(n,:)*2-1;
        nFy = elem(n,:)*2;
        Fn(nFx) = Fn(nFx)+Fe(1)/3;
        Fn(nFy) = Fn(nFy)+Fe(2)/3;
    end
end