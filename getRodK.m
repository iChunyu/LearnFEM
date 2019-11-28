% K = getRodKe(node,rod,E,A) calculates 2D stiffness matrix of rod elements
%     node --- coordinate of nodes in observer's frame
%     rod  --- node index of rod elements
%     E    --- elastic modulus of every element
%     A    --- area of every element

% More infomation:
%     The k-th row of 'node' represents the coordinate of k-th node
%     The k-th row of 'rod' represents the node of k-th rod

% XiaoCY 2019-11-27

%% main
function K = getRodK(node,rod,E,A)
    [Nnode,~] = size(node);
    [Nelem,~] = size(rod);
    K = zeros(2*Nnode);
    
    for k = 1:Nelem
        n1 = rod(k,1);
        n2 = rod(k,2);
        d1 = node(n1,:);
        d2 = node(n2,:);
        vec = d2-d1;
        L = sqrt(vec*vec');
        c = vec(1)/L;
        s = vec(2)/L;
        
        K(2*n1-1,2*n1-1) = K(2*n1-1,2*n1-1)+E(k)*A(k)/L*c^2;
        K(2*n1-1,2*n1) = K(2*n1-1,2*n1)+E(k)*A(k)/L*s*c;
        K(2*n1-1,2*n2-1) = K(2*n1-1,2*n2-1)-E(k)*A(k)/L*c^2;
        K(2*n1-1,2*n2) = K(2*n1-1,2*n2)-E(k)*A(k)/L*s*c;
        
        K(2*n1,2*n1-1) = K(2*n1,2*n1-1)+E(k)*A(k)/L*s*c;
        K(2*n1,2*n1) = K(2*n1,2*n1)+E(k)*A(k)/L*s^2;
        K(2*n1,2*n2-1) = K(2*n1,2*n2-1)-E(k)*A(k)/L*s*c;
        K(2*n1,2*n2) = K(2*n1,2*n2)-E(k)*A(k)/L*s^2;
        
        K(2*n2-1,2*n1-1) = K(2*n2-1,2*n1-1)-E(k)*A(k)/L*c^2;
        K(2*n2-1,2*n1) = K(2*n2-1,2*n1)-E(k)*A(k)/L*s*c;
        K(2*n2-1,2*n2-1) = K(2*n2-1,2*n2-1)+E(k)*A(k)/L*c^2;
        K(2*n2-1,2*n2) = K(2*n2-1,2*n2)+E(k)*A(k)/L*s*c;
        
        K(2*n2,2*n1-1) = K(2*n2,2*n1-1)-E(k)*A(k)/L*s*c;
        K(2*n2,2*n1) = K(2*n2,2*n1)-E(k)*A(k)/L*s^2;
        K(2*n2,2*n2-1) = K(2*n2,2*n2-1)+E(k)*A(k)/L*s*c;
        K(2*n2,2*n2) = K(2*n2,2*n2)+E(k)*A(k)/L*s^2;
    end
end