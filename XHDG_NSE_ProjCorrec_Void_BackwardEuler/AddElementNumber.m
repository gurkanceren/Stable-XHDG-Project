


function AddElementNumber(X,T)


for i=1:length(T)
    
    vertex1x=X(T(i,1),1);
    vertex1y=X(T(i,1),2);
    
    vertex2x=X(T(i,2),1);
    vertex2y=X(T(i,2),2);
    
    vertex3x=X(T(i,3),1);
    vertex3y=X(T(i,3),2);
    
    COG1=(vertex1x+vertex2x+vertex3x)/3;
    COG2=(vertex1y+vertex2y+vertex3y)/3;
    
    figure (1); hold on;
    
    %k=num2str(i);
    
    text(COG1,COG2,[num2str(i)]);
    
    hold on;
    
end
