% This function is to solve strain and stress
% [epsi,sigm] = solveStrAll(node,elem,A,d,E,mu,type)
%     epsi --- strain
%     sigm --- stress
%     node --- coordinate of nodes in observer's frame
%     elem --- node index of triangle elements
%     A    --- elements' area
%     d    --- nodes' displacement
%     E    --- elastic modulus of every element
%     mu   --- Poisson ratio
%     type --- type = 1 for plane stress problem
%              type = 2 for plane strain problem

% XiaoCY 2019-11-28

%% main
function [epsi,sigm] = solveStrAll(node,elem,A,d,E,mu,type)
    [Nnode,~] = size(node);
    [Nelem,~] = size(elem);
    strA = zeros(Nnode,7);
        
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
         
        b1 = y2-y3;
        b2 = y3-y1;
        b3 = y1-y2;
        c1 = x3-x2;
        c2 = x1-x3;
        c3 = x2-x1;
        
        BA = [ b1  0 b2  0 b3  0
            0 c1  0 c2  0 c3
            c1 b1 c2 b2 c3 b3]/2;
        
        idx = [ elem(n,1)*2-1
                elem(n,1)*2
                elem(n,2)*2-1
                elem(n,2)*2
                elem(n,3)*2-1
                elem(n,3)*2 ];
        dn = d(idx);
        
        epsiAn = BA*dn;
        sigmAn = D*epsiAn;
        strA(elem(n,:),1:3) = strA(elem(n,:),1:3)+repmat(epsiAn',3,1);
        strA(elem(n,:),4:6) = strA(elem(n,:),4:6)+repmat(sigmAn',3,1);
        strA(elem(n,:),7) = strA(elem(n,:),7)+A(n);
    end
    
    epsi = strA(:,1:3)./strA(:,7);
    sigm = strA(:,4:6)./strA(:,7);
end