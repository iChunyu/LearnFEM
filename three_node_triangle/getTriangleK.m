% K = getTriangleK(node,elem,E,mu,t,type) compute stiffness for plane problem
% Three-node triangle element
%     node --- coordinate of nodes in observer's frame
%     elem --- node index of triangle elements
%     E    --- elastic modulus of every element
%     mu   --- Poisson ratio
%     t    --- thickness
%     type --- type = 1 for plane stress problem
%              type = 2 for plane strain problem

% XiaoCY 2019-11-28

%% main
function [K,A] = getTriangleK(node,elem,E,mu,t,type)
    [Nnode,~] = size(node);
    [Nelem,~] = size(elem);
    K = zeros(2*Nnode);
    A = zeros(Nelem);
    
    switch type
        case 1
            D = [ 1   mu    0
                mu    1    0
                0    0  (1-mu)/2];
            D = E/(1-mu^2)*D;
        case 2
            E = E/(1-mu^2);
            mu = mu/(1-mu);
            D = [ 1   mu    0
                mu    1    0
                0    0  (1-mu)/2];
            D = E/(1-mu^2)*D;
        otherwise
            error('Wrong type flag!')
    end
    
    for n = 1:Nelem
        x1 = node(elem(n,1),1);
        y1 = node(elem(n,1),2);
        x2 = node(elem(n,2),1);
        y2 = node(elem(n,2),2);
        x3 = node(elem(n,3),1);
        y3 = node(elem(n,3),2);
        
%         a1 = x2*y3-x3*y2;
%         a2 = x3*y1-x1*y3;
%         a3 = x1*y2-x2*y1;
        An = det([1 x1 y1; 1 x2 y2; 1 x3 y3])/2;
        An = abs(An);
        A(n) = An;
        
        b1 = y2-y3;
        b2 = y3-y1;
        b3 = y1-y2;
        c1 = x3-x2;
        c2 = x1-x3;
        c3 = x2-x1;
        
        B = [ b1  0 b2  0 b3  0
            0 c1  0 c2  0 c3
            c1 b1 c2 b2 c3 b3]/An/2;
        Ke = B'*D*B*t*An;
        
        Index = zeros(1,6);
        Index([1 3 5]) = elem(n,:)*2-1;
        Index([2 4 6]) = elem(n,:)*2;
        K(Index,Index) = K(Index,Index)+Ke;
    end
end