function Vect_CutPtsEdge =  CountCutPoints(LSe,referenceElement)
%
% Vect_CutPtsEdge = CountCutPoints(LSe,referenceElement)
% The function only works for triangular elements. 
% Input: 
%     LSe: nodal values of the level-set function
%     referenceElement
% Output: 
%      Vect_CutPtsEdge: 3x1 vector. Each component is the number of times
%      that the interface cuts each edge in the triangle.


tol = 1e-20; 

nDeg = referenceElement.degree; 
Xe_ref = referenceElement.NodesCoord; 
Edges = referenceElement.faceNodes; 
[nEdges,nNodes] = size(Edges);

MaxLS_Edges = zeros(nEdges,1);

for k = 1:nEdges
    MaxLS_Edges(k) = max(abs(LSe(Edges(k,:))));
end

aux = find(MaxLS_Edges < tol);
if length(aux) == 1
    Edge0 = aux;
    Edge1 = Edge0+1; if Edge1 > nEdges, Edge1 = 1; end
    Edge2 = Edge0-1; if Edge2 == 0, Edge2 = nEdges; end
    
    P1 = Xe_ref(Edges(Edge1,1),:);
    P2 = Xe_ref(Edges(Edge1,end),:);
    t = linspace(0,1,3*p+1)';
    v = P2 - P1;
    pts = [P1(1) + t*v(1), P1(2) + t*v(2)];
    
    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,pts);
    N = N(:,:,1)'; 
    LS_aux = N*LSe;
    aux = LS_aux(1:end-1).*LS_aux(2:end);
    aux2 = find(abs(aux) < tol); aux(aux2) = zeros(size(aux2));
    aux_Edge1 = length(find(aux<=0));
    
    P1 = Xe_ref(Edges(Edge2,1),:);
    P2 = Xe_ref(Edges(Edge2,end),:);
    t = linspace(0,1)';
    v = P2 - P1;
    pts = [P1(1) + t*v(1), P1(2) + t*v(2)];
    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,pts);
    N = N(:,:,1)'; 
    LS_aux = N*LSe;
    aux = LS_aux(1:end-1).*LS_aux(2:end);
    aux2 = find(abs(aux) < tol); aux(aux2) = zeros(size(aux2));
    aux_Edge2 = length(find(aux<=0));
    
    if aux_Edge1 == 1 && aux_Edge2 == 1
        Vect_CutPtsEdge = [1,1,0];
    else
        Vect_CutPtsEdge = [1,1,1];
    end
    
elseif length(aux) == 2
    Vect_CutPtsEdge = [1,1,1];
else
    Vect_CutPtsEdge = zeros(nEdges,1);
    for i = 1:nEdges
        
        LS_aux = LSe(Edges(i,:));
        if all(LS_aux >= 0) || all(LS_aux <= 0)
            Vect_CutPtsEdge(i) = length(find(LS_aux == 0)) ; 
        else
            P1 = Xe_ref(Edges(i,1),:);
            P2 = Xe_ref(Edges(i,end),:);
            t = linspace(0,1)';
            v = P2 - P1;
            pts = [P1(1) + t*v(1), P1(2) + t*v(2)];
            
            N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,pts);
            N = N(:,:,1)'; 
            LS_aux = N*LSe;
            aux = find(abs(LS_aux) < 1e-16); 
            LS_aux(aux) = zeros(size(aux)); 
            aux = LS_aux(1:end-1).*LS_aux(2:end);
            aux = find(aux<=0);
            if length(aux) <= 1
                Vect_CutPtsEdge(i) = length(aux);
            else
                dif_aux = aux(2:end) - aux(1:end-1);
                ind_aux = find(abs(dif_aux) < 10);
                aux(ind_aux) = [];
                Vect_CutPtsEdge(i) = length(find(aux));
            end
        end
        
    end
    

    distVert = LSe(1:nEdges);
    aux = find(distVert == 0);
    for i = 1:length(aux)
        Edge1 = find(Edges(:,1) == aux(i));   LSe_1 = LSe(Edges(Edge1,:));
        Edge2 = find(Edges(:,end) == aux(i)); LSe_2 = LSe(Edges(Edge2,:));
       
        if LSe_1(2)*LSe_2(end-1) > 0 
            Vect_CutPtsEdge(Edge1) = Vect_CutPtsEdge(Edge1)-1;
            Vect_CutPtsEdge(Edge2) = Vect_CutPtsEdge(Edge2)-1;
        elseif Vect_CutPtsEdge(Edge1) < Vect_CutPtsEdge(Edge2)
            Vect_CutPtsEdge(Edge1) = Vect_CutPtsEdge(Edge1)-1;
        else
            Vect_CutPtsEdge(Edge2) = Vect_CutPtsEdge(Edge2)-1;
        end
    end
    aux = find(Vect_CutPtsEdge <0);
    Vect_CutPtsEdge(aux) = zeros(size(aux));
end





