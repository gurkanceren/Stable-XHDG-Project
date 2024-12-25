function [pc,qc]=computeElementsSolutionPressCorr(lambda_corr,UUc,QQc,Ufc,Qfc,T,F)
nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfFaceNodes = size(UUc{1},2)/3;

pc = zeros(nOfElements*nOfElementNodes,1);
qc = zeros(2*nOfElements*nOfElementNodes,1);

%Loop in elements
for ielem = 1:size(T,1)
    Fe = F(ielem,:);
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    pc((ielem-1)*nOfElementNodes+(1:nOfElementNodes)) = UUc{ielem}*lambda_corr(ind) + Ufc{ielem};
    qc((ielem-1)*2*nOfElementNodes+(1:2*nOfElementNodes)) = QQc{ielem}*lambda_corr(ind) + Qfc{ielem};
end
%
%disp('Hola');
