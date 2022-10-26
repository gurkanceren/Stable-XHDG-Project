function Elements = SetElements(T,LS,material,referenceElement)
%
% Elements = SetElements(T,LS,material,referenceElementl)
% 
% Input: 
%    T: matrix of element connectivities
%    LS: nodal level-set values
%    material: vector with two components, one for each subdomain. 
%              If material(i)=0, subdomain i is considered a void
%    referenceElement
% Output: 
%    Elements.D1: elements completely inside subdomain 1 (level set function > 0)
%    Elements.D2: elements completely inside subdomain 1 (level set function < 0)
%    Elements.Int: elements cut by the interface
%    Elements.Nodes: set of active nodes


tol = 1e-20; 

nDeg = referenceElement.degree;
Xe_ref = referenceElement.NodesCoord; 

nOfelemnodes=size(referenceElement.NodesCoord,1);

numel = size(T,1);


% In some cases, evaluating the level set function at the element's nodes
% is not enough. We will interpolate its value to evaluate it on a finer
% submesh, defined by points pts. 
npt = 5*nDeg+1;
x = linspace(-1,1,npt);
n = npt*(npt+1)/2;
pts = zeros(n,2);
ini = 0;
for i = 1:npt
    ind = ini + (1:npt-i+1);
    pts(ind,1) = x(1:npt-i+1)';
    pts(ind,2) = x(i);
    ini = i*(npt - 0.5*(i-1));
end
N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,pts);
N = N(:,:,1)'; 

ElemsD1  = zeros(numel,1); indD1 = 1;
ElemsD2  = zeros(numel,1); indD2 = 1;
ElemsInt = zeros(numel,1); indInt = 1;
Nodes = zeros(1,max(max(T)));
for ielem = 1:numel
    Te = T(ielem,1:nOfelemnodes);
    LSe = LS(Te);
    LS_aux = N*LSe; % Level set value at pts
    LS_aux(abs(LS_aux) < tol) = 0;
    LS_aux([1, npt,n]) = LSe(1:3); % the actual value of the level-set function is used at the vertices
    
    nodes1 = find(LS_aux >= 0); n1 = length(nodes1); 
    nodes2 = find(LS_aux <  0); n2 = length(nodes2); 
    if n1 < 3
        if max(abs(LS_aux(nodes1))) < 1e-15
            LS_aux(nodes1) = 0; 
        end
    end
    if n2 < 3
        if max(abs(LS_aux(nodes2))) < 1e-15
            LS_aux(nodes2) = 0; 
        end
    end
    
    if all(LS_aux >= 0)
        ElemsD1(indD1) = ielem;
        indD1 = indD1 + 1;
        if material(1) ~= 0
            Nodes(Te) = 1;
        end
    elseif  all(LS_aux <= 1e-12)
        ElemsD2(indD2) = ielem;
        indD2 = indD2 + 1;
        if material(2) ~= 0
            Nodes(Te) = 1;
        end
    else
        ElemsInt(indInt) = ielem;
        indInt = indInt + 1;
        Nodes(Te) = 1;
    end
end

ElemsD1  = ElemsD1(1:indD1-1);
ElemsD2  = ElemsD2(1:indD2-1);
ElemsInt = ElemsInt(1:indInt-1);
Nodes = find(Nodes == 1);

Elements.D1 = sort(ElemsD1);
Elements.D2 = sort(ElemsD2);
Elements.Int = sort(ElemsInt);
Elements.Nodes = Nodes;
