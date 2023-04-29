function SlopeTriangles (lhs,lerrors3,slopesFirstSegment)


% % errors = [0.017186	0.0030727	0.00044206
% %     0.0024447	0.00020942	1.8759e-05
% %     0.00031798	1.3653e-05	6.2886e-07
% % ];
% % 
% % hs = 0.5.^[3 4 5];
% % lerrors3 = log10(errors);
% % lhs = log10(hs);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%small triangle P3 
 
xcoord1=(lhs(3));
ycoord1=(lerrors3(3,3));
xcoord2=(lhs(2));
ycoord2=(lerrors3(3,3));
xcoord3=(lhs(2));
ycoord3=(lerrors3(2,3));


trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];


vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old


newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesFirstSegment(3)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)

%small triangle P2
xcoord1=(lhs(3));
ycoord1=(lerrors3(3,2));
xcoord2=(lhs(2));
ycoord2=(lerrors3(3,2));
xcoord3=(lhs(2));
ycoord3=(lerrors3(2,2));

trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];


vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old

newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesFirstSegment(2)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)

%small triangle P1
xcoord1=(lhs(3));
ycoord1=(lerrors3(3,1));
xcoord2=(lhs(2));
ycoord2=(lerrors3(3,1));
xcoord3=(lhs(2));
ycoord3=(lerrors3(2,1));

trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];

vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old

newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesFirstSegment(1)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slopesSecondSegment = (lerrors3(end-1,:)-lerrors3(end-2,:))/(lhs(end-1)-lhs(end-2));

slopesSecondSegment = ceil(slopesSecondSegment*100)/100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%small triangle P4 
 
xcoord1=(lhs(2));
ycoord1=(lerrors3(2,3));
xcoord2=(lhs(1));
ycoord2=(lerrors3(2,3));
xcoord3=(lhs(1));
ycoord3=(lerrors3(1,3));


trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];


vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old

newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesSecondSegment(3)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)

%small triangle P3
xcoord1=(lhs(2));
ycoord1=(lerrors3(2,2));
xcoord2=(lhs(1));
ycoord2=(lerrors3(2,2));
xcoord3=(lhs(1));
ycoord3=(lerrors3(1,2));

trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];


vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old

newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesSecondSegment(2)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)
%small triangle P2
xcoord1=(lhs(2));
ycoord1=(lerrors3(2,1));
xcoord2=(lhs(1));
ycoord2=(lerrors3(2,1));
xcoord3=(lhs(1));
ycoord3=(lerrors3(1,1));

trix=[xcoord1; xcoord2; xcoord3];
triy=[ycoord1; ycoord2; ycoord3];

vertices=[trix,triy].';
barycenter=mean(vertices,2);

newvertices=bsxfun(@plus, bsxfun(@minus,vertices,barycenter)*0.2, barycenter);
newbarycenter=mean(newvertices,2); %same as old

newvertices(:,4)=newvertices(:,1);
hold on 
plot(newvertices(1,:),newvertices(2,:),'k','LineWidth',2);

text(newvertices(1,2)+0.01,(newvertices(2,2)+0.1),num2str(slopesSecondSegment(1)),'HorizontalAlignment','left','FontSize',10,'LineWidth',1.5)
hold on
text((newvertices(1,1)+newvertices(1,2))/2,(newvertices(2,2))-0.2,num2str(1),'HorizontalAlignment','center','FontSize',10,'LineWidth',1.5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

