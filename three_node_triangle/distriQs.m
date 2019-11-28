% Fn = distriQs(node,inode,q) distribute concentrated force to nodes
%     node  --- coordinate of nodes in observer's frame
%     q     --- distributed force, contains [nodeIdx qx qy]

% XiaoCY 2019-11-28

%% main
function Fn = distriQs(node,q)
    [Nnode,~] = size(node);
    [Nq,~] = size(q);
    Fn = zeros(Nnode*2,1);
    
    for n = 1:Nq-1
        n1 = q(n,1);
        n2 = q(n+1,1);
        q1 = q(n,[2 3])';
        q2 = q(n+1,[2 3])';
        L = sqrt((node(n1,1)-node(n2,1))^2+(node(n1,2)-node(n2,2))^2);
        
        idx1 = [2*n1-1; 2*n1];
        idx2 = [2*n2-1; 2*n2];
        Fn(idx1) = Fn(idx1)+q1*L/3+q2*L/6;
        Fn(idx2) = Fn(idx2)+q1*L/6+q2*L/3;
    end
end