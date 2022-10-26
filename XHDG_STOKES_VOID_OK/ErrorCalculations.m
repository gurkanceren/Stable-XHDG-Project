
       
        %analytical solution
        u1_ex = @analyticalVelocityStokes_u1;
        u2_ex = @analyticalVelocityStokes_u2;

        
        %Relative error
        
        nOfD1Elements=length(Elements.D1);
        nOfD2Elements=length(Elements.D2);
        nOfCutElements=length(Elements.Int);
        vector=[Elements.D1];
        
        error_standard = zeros(nOfD1Elements,1);
        error_cut= zeros(nOfCutElements,1);
        error_standard_post=zeros(nOfD1Elements,1);
        error_cut_post=zeros(nOfCutElements,1);
        nOfElementNodes = size(referenceElement.NodesCoord,1);
        nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);


        
        for i = 1:nOfD1Elements          
            ielem=vector(i);
            ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
            error_standard(i) = computeL2Norm(referenceElement,X(T(ielem,:),:),(1:nOfElementNodes),u(ind,1),u1_ex).^2 ...
                +computeL2Norm(referenceElement,X(T(ielem,:),:),(1:nOfElementNodes),u(ind,2),u2_ex).^2;
                
            
        end
        
        Error_standard = sqrt(sum(error_standard));
        %disp(['XHDG Error on Standard Elements  = ', num2str(Error_standard)]);
        
        
        for i= 1:length(Elements.Int)
            ielem=Elements.Int(i);
            ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
            error_cut(i) = computeL2Norm_cut(LS(T(ielem,:)),referenceElement,X(T(ielem,:),:),(1:nOfElementNodes),u(ind,1),u1_ex).^2 ...
                +computeL2Norm_cut(LS(T(ielem,:)),referenceElement,X(T(ielem,:),:),(1:nOfElementNodes),u(ind,2),u2_ex).^2 ;
            
        end
        Error_cut = sqrt(sum(error_cut));
        %disp(['XHDG Error over Cut Elements = ', num2str(Error_cut)]);
                
        Error_total=sqrt(Error_cut^2+Error_standard^2);
        disp(['XHDG Error Overall Domain = ', num2str(Error_total)]);

              
        
        for i = 1:nOfD1Elements         
            ielem=vector(i);
            ind = (ielem-1)*nOfElementNodes_star+1:ielem*nOfElementNodes_star;
            Xold=X(T(ielem,:),:);
            Xe=shapeFunctionss*Xold;
            error_standard_post(i) = computeL2Norm(referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,1),u1_ex).^2 ...
                +computeL2Norm(referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,2),u2_ex).^2;
                       
        end

         
        Error_standard_post = sqrt(sum(error_standard_post));
        %disp(['XHDG Error on Standard Elements PostProcessed  = ', num2str(Error_standard_post)]);
        
        
        for i= 1:length(Elements.Int)
            ielem=Elements.Int(i);
            ind = (ielem-1)*nOfElementNodes_star+1:ielem*nOfElementNodes_star;
            Xold=X(T(ielem,:),:);
            Xe=shapeFunctionss*Xold;
            LSe_star = EvaluateLS(Xe,example);            
            error_cut_post(i) = computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,1),u1_ex).^2 ... 
                + computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,2),u2_ex).^2;
            
        end

        Error_cut_post = sqrt(sum(error_cut_post));
        %disp(['XHDG Error over Cut Elements PostProcessed= ', num2str(Error_cut_post)]);
                
        Error_total_post=sqrt(Error_cut_post^2+Error_standard_post^2);
        disp(['XHDG Error Overall Domain Postprocessed = ', num2str(Error_total_post)]);
         

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        