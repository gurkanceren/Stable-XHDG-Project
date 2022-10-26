function lap=analyticaLaplacianBimaterialCircleExample(x)

n=size(x,1);

lap = zeros(n,1);

R = 0.4;
a = 0.4; mu1=1; mu2=5;
A = a*(2*mu1+mu2)/(2*mu2);
B = a^3*(mu2-2*mu1)/(2*mu2);

for i=1:n
  r = norm(x(i,:));
  if r>R
      dudr = 2*A*r-2*B/r^3; %derivative of u with respect to r
      d2udr2 = 2*A+6*B/r^4; %2nd derivative of u with respect to r
  else
      dudr = 2*r;
      d2udr2 = 2;
  end
  
  drdx = x(i,1)/r; drdy= x(i,2)/r; %derivative of r with respect to x and y
  d2rdx2 = -(x(i,1)^2)/r^3 + 1/r; %2nd derivative of r with respect to x and y
  d2rdy2 = -(x(i,2)^2)/r^3 + 1/r;
  
  lap(i) = d2udr2*drdx^2+dudr*d2rdx2 + d2udr2*drdy^2+dudr*d2rdy2;
  
end
    