function Pts = FindPointsInterface(LSe, referenceElement)
%
% Pts = FindPointsInterface(LSe, referenceElement)
%
% Given a referenceElement and the value of the level set function at the
% nodes, this function computes a set of p+1 points (p is the interpolation
% degree) that define a polynomial approximating the interface. 
% This function works for elements that are cut by the interface in an
% "standard way". 

Xe_aux = [-1,-1; 1,-1; -1,1];
LS_vertices = LSe(1:3);
aux = find(LS_vertices == 0);
if length(aux) == 2
    P1 = Xe_aux(aux(1),:);
    P2 = Xe_aux(aux(2),:);
else
    aux = find(LS_vertices > 0); n1 = length(aux);
    if n1 == 1
        node1 = find(LS_vertices > 0);
    else
        node1 = find(LS_vertices <= 0);
    end
    node2 = node1+1; if node2 == 4, node2 = 1; end
    node3 = node1-1; if node3 == 0, node3 = 3; end
    
    P1_aux = Xe_aux(node1,:); d1_aux = LS_vertices(node1);
    P2_aux = Xe_aux(node2,:); d2_aux = LS_vertices(node2);
    P1 = FindPointLS(0,LSe,P1_aux,d1_aux,P2_aux,d2_aux,referenceElement);
    if any(P1 < min(min(Xe_aux))) || any(P1 > max(max(Xe_aux)))
        PM = (P1_aux + P2_aux)/2;
        N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,PM);
        N = N(:,:,1)'; 
        dM = N*LSe;
        if d1_aux*dM < 0
            P1 = FindPointLS(0,LSe,P1_aux,d1_aux,PM,dM,referenceElement);
        else
            P1 = FindPointLS(0,LSe,PM,dM,P2_aux,d2_aux,referenceElement);
        end
    end
    
    P1_aux = Xe_aux(node1,:); d1_aux = LS_vertices(node1);
    P2_aux = Xe_aux(node3,:); d2_aux = LS_vertices(node3);
    P2 = FindPointLS(0,LSe,P1_aux,d1_aux,P2_aux,d2_aux,referenceElement);
    if any(P2 < min(min(Xe_aux))) || any(P2 > max(max(Xe_aux)))
        PM = (P1_aux + P2_aux)/2;
        N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,PM);
        N = N(:,:,1)'; 
        dM = N*LSe;
        if d1_aux*dM < 0
            P2 = FindPointLS(0,LSe,P1_aux,d1_aux,PM,dM,referenceElement);
        else
            P2 = FindPointLS(0,LSe,PM,dM,P2_aux,d2_aux,referenceElement);
        end
    end
end


Pts = FindPointspP(P1,P2,LSe,referenceElement);





%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------



function res = FindPointspP(P1,P2,LSe,referenceElement)
%
% res = FindPointspP(P1,P2,elem,nen,p,LSe,Xe_ref,FShapeFunc)


Xe_ref = referenceElement.NodesCoord; 
p = referenceElement.degree;

if p == 1
    res = [P1; P2];
else
    v = P2-P1;
    n = [-v(2), v(1)];

    y1 = min(Xe_ref(:,2));
    aux = find(abs(Xe_ref(:,2) - y1) < 1e-6);
    aux = sort(Xe_ref(aux,1));
    x1 = min(aux); x2 = max(aux);
    lambda = aux/(x2-x1) - x1/(x2-x1);

    res = zeros(p+1,2);
    res(1,:) = P1;
    res(p+1,:) = P2;
    for i = 2:p
        P_aux = P1  + lambda(i)*v;
        
        P0 = FindPointLS(0,LSe,P_aux-n/2,[],P_aux+n/2,[],referenceElement);
        res(i,:) = P0;
    end
end


