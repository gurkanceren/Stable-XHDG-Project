
function [Nzg,N_enrc,Iprime,Pphy,g]=modifiedReferenceElementInterface(referenceElement,PtsInt,Xe)

    p=referenceElement.degree;
    g=length(referenceElement.IPweights1d);
    
    PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
    
    zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
    Nzg = shapeFunctions(:,:,1)';  % 2d shape functions on the refence element at integration points
    
    Nzgl=Nzg*1;
    Nzgr=Nzg*-1;
    
    NzgH=(Nzgl+Nzgr)/2;
    Nzg=[Nzg  NzgH];
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
    NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
    
    NPtsIntl=NPtsInt*1;
    NPtsIntr=NPtsInt*-1;
    
    NPtsIntH=(NPtsIntl+NPtsIntr)/2;
    NPtsInt=[NPtsInt NPtsIntH];
    
    Pphy=NPtsInt*Xe; %p+1 interface nodes on the REAL element
       
    Iprime=referenceElement.N1d*Pphy;   %*Pphy;





































% %     p=referenceElement.degree;
% %     g=length(referenceElement.IPweights1d);
% %   
% %     
% %     
% %     PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
% %     
% %     zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
% %     
% %     shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
% %     Nzg = shapeFunctions(:,:,1)';  % 2d shape functions on the refence element at integration points
% %     
% %     shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
% %     NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
% %     
% %     Pphy=NPtsInt*Xe; %p+1 interface nodes on the REAL element
% %        
% %     Iprime=referenceElement.N1dxi*Pphy;   %*Pphy;