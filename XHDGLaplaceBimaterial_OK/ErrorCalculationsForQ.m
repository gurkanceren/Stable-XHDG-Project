

        %q0 = @analiticalSolutionq;
        
        %Relative error
        
        nOfD1Elements=length(Elements.D1);
        nOfD2Elements=length(Elements.D2);
        nOfCutElements=length(Elements.Int);
        vector=[Elements.D1;Elements.D2];
        
        error_standard = zeros(nOfD1Elements+nOfD2Elements,1);
        error_cut= zeros(nOfCutElements,1);
        
        for i = 1:nOfD1Elements+nOfD2Elements          
            ielem=vector(i);
            ind = (ielem-1)*(4*nOfElementNodes) + (1:4*nOfElementNodes);
            error_standard_q(i) = computeL2Normq(referenceElement,X(T(ielem,:),:),(1:nOfElementNodes),q(ind(1:2*nOfElementNodes)));
            
        end
        
        Error_standard_q = sqrt(sum(error_standard_q.^2));
        disp(['XHDG Error on Standard Elements (Q)  = ', num2str(Error_standard_q)]);
        
        
        for i= 1:length(Elements.Int)
            ielem=Elements.Int(i);
            ind = (ielem-1)*(4*nOfElementNodes) + (1:4*nOfElementNodes);
            error_cut_q(i) = computeL2Norm_cutq(LS(T(ielem,:)),referenceElement,X(T(ielem,:),:),(1:4*nOfElementNodes),q(ind));
            
        end
        Error_cut_q = sqrt(sum(error_cut_q.^2));
        disp(['XHDG Error over Cut Elements (Q) = ', num2str(Error_cut_q)]);
                
        Error_total_q=sqrt(Error_cut_q^2+Error_standard_q^2);
        disp(['XHDG Error Overall Domain (Q) = ', num2str(Error_total_q)]);