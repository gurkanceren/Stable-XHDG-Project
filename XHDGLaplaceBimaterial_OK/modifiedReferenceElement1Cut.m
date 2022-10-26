
function [Nx_cutf_old,Nx_cutf,N_cutf,wgp_f,glim]=modifiedReferenceElement1Cut(referenceElement,LSe,nodes)


[zgp_f,wgp_f,n1_f,n2_f,IntPt] = ModifyQuadrature1D(LSe(nodes,1),referenceElement);

shapeFunctions_f=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord1d,zgp_f);

N_cutf=shapeFunctions_f(:,:,1)';
Nl_f=N_cutf(1:n1_f,:)*-1;
Nr_f=N_cutf(n1_f+1:n2_f+n1_f,:)*1;
NH_f=[Nl_f;Nr_f];
N_cutf=[N_cutf,NH_f];



Nx_cutf_old=shapeFunctions_f(:,:,2)';
Nxl_f=Nx_cutf_old(1:n1_f,:)*-1;
Nxr_f=Nx_cutf_old(n1_f+1:n2_f+n1_f,:)*1;
NHx_f=[Nxl_f;Nxr_f];
Nx_cutf=[Nx_cutf_old,NHx_f];


glim=n1_f+n2_f;
