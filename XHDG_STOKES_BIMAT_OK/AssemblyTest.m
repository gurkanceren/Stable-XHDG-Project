clear all, close all, %home
restoredefaultpath, setpath
    
%% Load computational mesh 
meshName = 'mesh6_P2.dcm';
if all(meshName(end-2:end)=='dcm'),GenerateMatFileFromEZ4U(['Meshes/' meshName]); end
load(['Meshes/' meshName(1:end-3) 'mat']); 
nOfElements = size(T,1); nOfElementNodes = size(T,2);
figure(1),clf,plotMesh(X,T)

referenceElement = createReferenceElement(1,elemInfo.nOfNodes);

[numel,nen] = size(T);
npt = size(X,1); 

Ke = rand(nen); 


tic
K1 = spalloc(npt,npt,nen^2*npt); 
for ielem = 1:numel
    Te = T(ielem,:);
    K1(Te,Te) = K1(Te,Te) + Ke; 
end
toc



tic
ind = 1; n = nen^2*numel;
ind_i  = zeros(1,n); ind_j  = zeros(1,n);
coef_K = zeros(1,n);
for ielem = 1:numel
    Te = T(ielem,:);
    for irow = 1:nen
        for icol = 1:nen
            ind_i(ind)  = Te(irow);
            ind_j(ind)  = Te(icol);
            coef_K(ind) = Ke(irow,icol);
            ind = ind+1;
        end
    end
end
ind_i  = ind_i(1:ind-1);
ind_j  = ind_j(1:ind-1);
coef_K = coef_K(1:ind-1);
K2 = sparse(ind_i,ind_j,coef_K);
toc

tic
ind = 1:nen;
n = nen^2*numel;
ind_i  = zeros(1,n); ind_j  = zeros(1,n);
coef_K = zeros(1,n);
for ielem = 1:numel
    Te = T(ielem,:);
    for irow = 1:nen
        ind_i(ind) = Te(irow);
        ind_j(ind) = Te';
        coef_K(ind) = Ke(irow,:);
        ind = ind + nen;
    end
end
K3 = sparse(ind_i,ind_j,coef_K);
toc


tic
n = nen^2*numel;
ind_i  = zeros(1,n); ind_j  = zeros(1,n); coef_K = zeros(1,n);
for ielem = 1:numel
    Te = T(ielem,:);
    vect = (ielem-1)*nen^2 + (1:nen^2); 
    ind_i( vect ) = repmat(Te',nen,1); 
    ind_j( vect ) = reshape(repmat(Te,nen,1),nen^2,1); 
    coef_K( vect) = reshape(Ke,nen^2,1); 
end
K4 = sparse(ind_i,ind_j,coef_K);
toc


tic
n = nen^2*numel;
ind_i  = zeros(1,n);
ind_j  = zeros(1,n);
coef_K = zeros(1,n);
for ielem = 1:numel
    Te = T(ielem,:);
    
    [rows,cols] = meshgrid(Te,Te); 
    rows = reshape(rows,nen^2,1); 
    cols = reshape(cols,nen^2,1); 
    Ke_aux = reshape(Ke',nen^2,1); 
    
    vect = (ielem-1)*nen^2 + (1:nen^2); 
    ind_i( vect ) = rows; 
    ind_j( vect ) = cols; 
    coef_K( vect) = Ke_aux; 
end
K5 = sparse(ind_i,ind_j,coef_K);
toc

dif_12 = max(max(abs(K1 - K2)))
dif_23 = max(max(abs(K1 - K3)))
dif_14 = max(max(abs(K1 - K4)))
dif_15 = max(max(abs(K1 - K5)))
