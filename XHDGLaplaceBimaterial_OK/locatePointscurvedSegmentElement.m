function newPoints = locatePointscurvedSegmentElement(coord1D,Points)
%newPoints = locatePoints1DcurvedElement(coord1D,Points)
%Given a curved segment defined by p+1 point in Points ( (p+1)x nsd )
% and the coordinates in the 1D reference element
%returns the coordinates of the p+1 points on the same segment with
%the distances (in arclength) given by coord1D.
%The 1st and last point are assumed to coorrespond to the end points of the
%segment

%The same reference element is assumed for the data and the output. It the
%data corresponds to equally spaced nodes, better mofidy this variable
coord1D0 = coord1D;

p=length(coord1D0)-1; %degree

%Computation of the lengths
l = zeros(p+1,1);
for i=2:p+1
   [z,w]=gaussLegendre(p,-1,coord1D0(i)); %quadrature for [-1,s_i]
   shapeFunctions = computeShapeFunctionsAtPoints(p,coord1D0,z);
   Nxi = shapeFunctions(:,:,2)';
   sum = 0;
   for k=1:length(w)
       Nxi_k = Nxi(k,:); %functions derivatives at k-th integration point
       velocity = Nxi_k*Points;
       veloNorm = norm(velocity);
       sum = sum +  veloNorm*w(k);
   end
   l(i) = sum;   
end
%Normalization for [-1,1]
l = 2*l/l(end)-1;

%Computation of new coordinates
%Basis functions for points l evaluated at coord1D
N = computeOnlyShapeFunctionsAtPoints1D(coord1D,p,l); 
newPoints = N*Points;
%The 1st and last point are assumed to coorrespond to the end points
newPoints(1,:)=Points(1,:); newPoints(end,:)=Points(end,:);



