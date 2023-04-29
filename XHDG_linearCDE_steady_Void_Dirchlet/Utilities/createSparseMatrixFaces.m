function A=createSparseMatrixFaces(nDOFface,intFaces,extFaces,F,c)
% nDOFface: number of degrees of freedom for each face
% matrix with c value for all non-null coefficients

nOfInteriorFaces=size(intFaces,1);
nOfExteriorFaces=size(extFaces,1);
i = []; j=[];
for f=1:nOfInteriorFaces;
   ii=f*ones(1,5);
   ielem1=intFaces(f,1);
   ielem2=intFaces(f,3);
   jj=union(F(ielem1,:),F(ielem2,:));
   i = [i ii]; j=[j jj];
end
for f=1:nOfExteriorFaces
   ii=f*ones(1,3)+nOfInteriorFaces;
   ielem1=extFaces(f,1);
   jj=F(ielem1,:);
   i = [i ii]; j=[j jj];
end

n = nOfInteriorFaces+nOfExteriorFaces;
A = sparse(i,j,(c+1.e-300)*ones(1,length(i)),n,n);

%Convert to block sparse matrix, with blocks nDOFface^2
A = repmat(A,nDOFface,nDOFface);
p = [1:n:n*nDOFface]'; 
p = bsxfun(@plus,repmat(p,1,n),[0:n-1]); 
p = reshape(p,n*nDOFface,1);
A = A(p,p);




