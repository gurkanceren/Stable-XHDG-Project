
 
nOfNodes=size(X,1);
nOfElements=size(T,1);
nOfElementNodes=size(T,2);
x_to_search=zeros(nOfElements*3,1);
x=X(:,1);
for i=1:nOfElements
telem=T(i,:);
xelem=x(T(i,:));
ind = (i-1)*(3) + (1:3);
x_to_search(ind,1)=xelem(1:3);
end
b=0.2;
[m,i] = min(abs(x_to_search-b));
a = x_to_search(i);

i = find(x<=a);
j = find(x>a);

x(i) = ((b+1)/(a+1))*(x(i)+1)-1;
x(j) = ((b-1)/(a-1))*(x(j)-1)+1;

X(:,1)=x;
