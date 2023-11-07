function XeNew = optimalPointsIsoparametric2D(Xe,p,XeRef)
%
% XeNew = optimalPointsIsoparametric2D(Xe,p,XeRef)
% 
% Xe: nodal coordinates of the element (x-y)
% p: degree of interpolation
% XeRef: nodal coordinates of the element (xi-eta)
%
% NOTE. Reference triangle is [0,0; 1,0; 0,1];
%

% 1D nodal distribution (at eta=0)
sc = XeRef(XeRef(:,2)==0,1);
sc = sort(sc);
sc1 = sort(sc,'descend');

% Boundary nodes on local coordiantes (xi-eta)
sb = [sc            , 0*sc; 
      sc1(2:end)    , sc(2:end); 
      0*sc1(2:end-1), sc1(2:end-1)];
  
% Boundary nodes on cartesian coordinates (x-y)
xp = 0*sb;
for i=1:size(sb,1)
    for j=1:size(XeRef,1)
        if norm(sb(i,:)-XeRef(j,1:2)) < 1.e-8
            xp(i,:) = Xe(j,:); 
        end
    end
end

V = zeros(size(sb,1),size(sb,1));
for i=1:size(sb,1)
    V(:,i) = optimalPointsIsoparametricBlendingBoundary2D(sb,i,p,sc); 
end
alpha = V\xp;

V = zeros(size(XeRef,1),size(sb,1));
for i=1:size(sb,1)
    V(:,i) = optimalPointsIsoparametricBlendingBoundary2D(XeRef(:,1:2),i,p,sc); 
end
XeNew = V*alpha;


%______________________________________________________________
function W = optimalPointsIsoparametricBlendingBoundary2D(Z,m,p,s)
%
% W = optimalPointsIsoparametricBlendingBoundary2D(Z,m,p,s)
%

if m == 1
    % First vertex ([0,0] in xi-eta)
    W = 1-Z(:,1)-Z(:,2);
elseif m == p+1
    % Second vertex ([1,0] in xi-eta)
    W = Z(:,1);
elseif m == 2*p+1
    % Third vertex ([0,0] in xi-eta)
    W = Z(:,2);
else
    % Edge nodes
    C = ones(p-1,p-1);
    C(:,1) = s(2:end-1).*(1-s(2:end-1));
    for i=2:p-1
        C(:,i) = C(:,i-1).*s(2:end-1); 
    end
    C = inv(C);
    if m < p+1 
        % First edge 
        ZC = Z; 
        ind = m-1;
    elseif m < 2*p+1 
        % Second edge 
        ZC = [Z(:,2),1-Z(:,1)-Z(:,2)]; 
        ind = m-p-1;
    else 
        % Third edge 
        ZC = [1-Z(:,1)-Z(:,2),Z(:,1)]; 
        ind = m-2*p-1; 
    end
    W = Z(:,1)*0;
    for i=1:p-1
        W = W + C(i,ind)*ZC(:,1).^i; 
    end
    W = W.*(1-ZC(:,1)-ZC(:,2));
end


