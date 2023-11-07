


if ~isempty(d1) %d1 Elements
    
    PtsInt=[];
    FaceInfo=[];
    
    [Auq,Alu,Alq,Auu,Aqq,All,fe]=xhdg_matrix_d1(referenceElement,mu1,PtsInt,FaceInfo,tau,Xe,LSe,iElem);
    
    % Elemental mapping
    Aqu=-Auq';
    Aul = -Alu'; Aql = Alq';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    fUQ= A\[fe;zeros(2*nOfElementNodes,1)];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); 
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); 
    
elseif ~isempty(d2) %d2 Elements
    
    PtsInt=[];
    FaceInfo=[];
    
    [Auq,Alu,Alq,Auu,Aqq,All,fe]=xhdg_matrix_d2 (referenceElement,mu2,PtsInt,FaceInfo,tau,Xe,LSe,iElem);
    
    % Elemental mapping
    Aqu=-Auq';
    Aul = -Alu'; Aql = Alq';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    fUQ= A\[fe;zeros(2*nOfElementNodes,1)];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); 
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end);