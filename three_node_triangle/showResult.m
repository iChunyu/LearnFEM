% showResult(node,elem,rst) plots result in color map
%     node --- coordinate of nodes in observer's frame
%     elem --- node index of triangle elements
%     rst  --- result to be showed

% XiaoCY 2019-11-28

%% main
function showResult(node,elem,rst)
    [Nelem,~] = size(elem);
    x = zeros(3,Nelem);
    y = zeros(3,Nelem);
    c = zeros(3,Nelem);
    
    figure
    for n = 1:Nelem
        xn = node(elem(n,:),1);
        yn = node(elem(n,:),2);
        cn = rst(elem(n,:));
        
        x(:,n) = xn;
        y(:,n) = yn;
        c(:,n) = cn;
    end
    patch(x,y,c,'Facecolor','interp','EdgeColor','interp')
    grid on
    colorbar
    axis equal
end