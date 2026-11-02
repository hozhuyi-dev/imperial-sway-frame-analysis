classdef FRAME
    
    % Ho Zhu Yi
    % 01082991
    % 28 February 2017
    
    % All dimensions in N and mm.
    
    properties
        % Constants
        e
        r
        ea
        ei
        volume
        
        % Nodes coordinates and it's dofs
        ncoords
        ndofs
        ncoordsnew
        ncoordsplot
        
        % Nodes connected to each element
        enodes
        
        % Free and restrained dofs
        dfree
        drestrained
        
        % Number of nodes and elements
        nodes
        elements
        
        % Stiffness matrix, force, deflection
        k
        f
        u
        
        % Maximum lateral sway
        uxmax
        
        % Equilibrium checks
        exthor
        extver
        check
        
        % Buckling checks: Euler buckling, axial load, margin of safety
        pcr
        axiload
        mos
        mosmin
    end
    
    methods
        % Creation of a variable with this class
        function obj=FRAME(E,R,NCOORDS,ENODES,DFREE,DRESTRAINED)
            obj.e=E;
            obj.r=R;
            obj.ncoords=NCOORDS;
            obj.enodes=ENODES;
            obj.dfree=DFREE;
            obj.drestrained=DRESTRAINED;
            
            % Calculate ea and ei
            obj.ea=obj.e.*(pi*(obj.r.^2));
            obj.ei=obj.e.*(pi*(obj.r.^4)/4);
            
            % Calculate the number of nodes and elements
            obj.nodes=size(obj.ncoords,1);
            obj.elements=size(obj.enodes,1);
            
            % Defining the numbering of dofs
            obj.ndofs=zeros(obj.nodes,3);
            for i=1:obj.nodes
                obj.ndofs(i,:)=[3*i-2,3*i-1,3*i];
            end
        end
        
        % Calculate the total volume
        function obj=Volume(obj)
            obj.volume=0;
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                Le=sqrt((x2-x1)^2+(y2-y1)^2 );
                obj.volume=obj.volume+(Le*pi*(obj.r(i)^2));
            end
        end
        
        % Solve using specific u restrained and f
        function obj=Solve(obj,U,F)
            % Construct stiffness matrix k
            obj.k=zeros(3*(obj.nodes));
            for i=1:obj.elements % loop through all elements & build stiffness matrix
                n1=obj.enodes(i,1); n2=obj.enodes(i,2); % identify element node numbers
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2); % element node 1 - x,y coordinates
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2); % element node 2 - x,y coordinates
                dof11=obj.ndofs(n1,1); dof12=obj.ndofs(n1,2); dof13=obj.ndofs(n1,3); % element node 1 - dofs
                dof21=obj.ndofs(n2,1); dof22=obj.ndofs(n2,2); dof23=obj.ndofs(n2,3); % element node 2 - dofs
                alpha=atan2(y2-y1,x2-x1); % angle of inclination of element to POSITIVE direction of the horizontal x axis
                c=cos(alpha); c2=c*c; s=sin(alpha); s2=s*s; cs=c*s; % angle parameters
                Le=sqrt((x2-x1)^2+(y2-y1)^2); % element length
                ka=obj.ea(i)/Le; % generic stretching stiffness term
                kb=obj.ei(i)/(Le*Le*Le); % generic bending stiffness term
                
                % Updating global stiffness matrix [K] coefficients
                % Note that for each element you still have to 'think locally'
                % Row 1 - element dof11 - u1 of element
                obj.k(dof11,dof11)=obj.k(dof11,dof11)+ka*c2+12*kb*s2;      % Col 1 - element dof11 - u1 of element
                obj.k(dof11,dof12)=obj.k(dof11,dof12)+ka*cs-12*kb*cs;      % Col 2 - element dof12 - v1 of element
                obj.k(dof11,dof13)=obj.k(dof11,dof13)-6*Le*kb*s;             % Col 3 - element dof13 - th1 of element
                obj.k(dof11,dof21)=obj.k(dof11,dof21)-ka*c2-12*kb*s2;      % Col 4 - element dof21 - u2 of element
                obj.k(dof11,dof22)=obj.k(dof11,dof22)-ka*cs+12*kb*cs;      % Col 5 - element dof22 - v2 of element
                obj.k(dof11,dof23)=obj.k(dof11,dof23)-6*Le*kb*s;             % Col 6 - element dof23 - th2 of element
                
                % Row 2 - element dof12 - v1 of element
                obj.k(dof12,dof11)=obj.k(dof12,dof11)+ka*cs-12*kb*cs;      % Col 1 - element dof11 - u1 of element
                obj.k(dof12,dof12)=obj.k(dof12,dof12)+ka*s2+12*kb*c2;      % Col 2 - element dof12 - v1 of element
                obj.k(dof12,dof13)=obj.k(dof12,dof13)+6*Le*kb*c;             % Col 3 - element dof13 - th1 of element
                obj.k(dof12,dof21)=obj.k(dof12,dof21)-ka*cs+12*kb*cs;      % Col 4 - element dof21 - u2 of element
                obj.k(dof12,dof22)=obj.k(dof12,dof22)-ka*s2-12*kb*c2;      % Col 5 - element dof22 - v2 of element
                obj.k(dof12,dof23)=obj.k(dof12,dof23)+6*Le*kb*c;             % Col 6 - element dof23 - th2 of element
                
                % Row 3 - element dof13 - th1 of element
                obj.k(dof13,dof11)=obj.k(dof13,dof11)-6*Le*kb*s;             % Col 1 - element dof11 - u1 of element
                obj.k(dof13,dof12)=obj.k(dof13,dof12)+6*Le*kb*c;             % Col 2 - element dof12 - v1 of element
                obj.k(dof13,dof13)=obj.k(dof13,dof13)+4*Le*Le*kb;            % Col 3 - element dof13 - th1 of element
                obj.k(dof13,dof21)=obj.k(dof13,dof21)+6*Le*kb*s;             % Col 4 - element dof21 - u2 of element
                obj.k(dof13,dof22)=obj.k(dof13,dof22)-6*Le*kb*c;             % Col 5 - element dof22 - v2 of element
                obj.k(dof13,dof23)=obj.k(dof13,dof23)+2*Le*Le*kb;            % Col 6 - element dof23 - th2 of element
                
                % Row 4 - element dof21 - u2 of element
                obj.k(dof21,dof11)=obj.k(dof21,dof11)-ka*c2-12*kb*s2;      % Col 1 - element dof11 - u1 of element
                obj.k(dof21,dof12)=obj.k(dof21,dof12)-ka*cs+12*kb*cs;      % Col 2 - element dof12 - v1 of element
                obj.k(dof21,dof13)=obj.k(dof21,dof13)+6*Le*kb*s;             % Col 3 - element dof13 - th1 of element
                obj.k(dof21,dof21)=obj.k(dof21,dof21)+ka*c2+12*kb*s2;      % Col 4 - element dof21 - u2 of element
                obj.k(dof21,dof22)=obj.k(dof21,dof22)+ka*cs-12*kb*cs;      % Col 5 - element dof22 - v2 of element
                obj.k(dof21,dof23)=obj.k(dof21,dof23)+6*Le*kb*s;             % Col 6 - element dof23 - th2 of element
                
                % Row 5 - element dof22 - v2 of element
                obj.k(dof22,dof11)=obj.k(dof22,dof11)-ka*cs+12*kb*cs;      % Col 1 - element dof11 - u1 of element
                obj.k(dof22,dof12)=obj.k(dof22,dof12)-ka*s2-12*kb*c2;      % Col 2 - element dof12 - v1 of element
                obj.k(dof22,dof13)=obj.k(dof22,dof13)-6*Le*kb*c;             % Col 3 - element dof13 - th1 of element
                obj.k(dof22,dof21)=obj.k(dof22,dof21)+ka*cs-12*kb*cs;      % Col 4 - element dof21 - u2 of element
                obj.k(dof22,dof22)=obj.k(dof22,dof22)+ka*s2+12*kb*c2;      % Col 5 - element dof22 - v2 of element
                obj.k(dof22,dof23)=obj.k(dof22,dof23)-6*Le*kb*c;             % Col 6 - element dof23 - th2 of element
                
                % Row 6 - element dof23 - th2 of element
                obj.k(dof23,dof11)=obj.k(dof23,dof11)-6*Le*kb*s;             % Col 1 - element dof11 - u1 of element
                obj.k(dof23,dof12)=obj.k(dof23,dof12)+6*Le*kb*c;             % Col 2 - element dof12 - v1 of element
                obj.k(dof23,dof13)=obj.k(dof23,dof13)+2*Le*Le*kb;            % Col 3 - element dof13 - th1 of element
                obj.k(dof23,dof21)=obj.k(dof23,dof21)+6*Le*kb*s;             % Col 4 - element dof21 - u2 of element
                obj.k(dof23,dof22)=obj.k(dof23,dof22)-6*Le*kb*c;             % Col 5 - element dof22 - v2 of element
                obj.k(dof23,dof23)=obj.k(dof23,dof23)+4*Le*Le*kb;            % Col 6 - element dof23 - th2 of element
            end
            
            % Initiate U and specify restrained dof
            obj.u=zeros(3*obj.nodes,1);
            obj.u(obj.drestrained)=U;
            
            % Construct the global nodal force vector and calculate the
            % total external horizontal and vertical force to be used in
            % equilibrium checks
            obj.f=zeros(3*obj.nodes,1);
            obj.exthor=0;
            obj.extver=0;
            for i=1:size(F,1)
                obj.f(F(i,1))=F(i,2);
                if rem(F(i,1),3)==1
                    obj.exthor=obj.exthor+F(i,2);
                elseif rem(F(i,1),3)==2
                    obj.extver=obj.extver+F(i,2);
                end
            end
            
            % Solve for the unknown nodal dofs
            obj.u(obj.dfree)=obj.k(obj.dfree,obj.dfree)\(obj.f(obj.dfree)-obj.k(obj.dfree,obj.drestrained)*obj.u(obj.drestrained));
            
            % Solve for the unknown reactions
            obj.f(obj.drestrained)=obj.k(obj.drestrained,obj.dfree)*obj.u(obj.dfree)+obj.k(obj.drestrained,obj.drestrained)*obj.u(obj.drestrained);
            
            % Updating nodal coordinates after deformation
            obj.ncoordsnew=zeros(size(obj.ncoords));
            for i=1:obj.nodes
                obj.ncoordsnew(i,1)=obj.ncoords(i,1)+obj.u(3*i-2);
                obj.ncoordsnew(i,2)=obj.ncoords(i,2)+obj.u(3*i-1);
            end
        end
        
        % Find the maximum lateral sway
        function obj=MaximumLateralSway(obj)
            ux=obj.u(1:3:end);
            obj.uxmax=max(ux);
        end
        
        % Plot the final structure
        function obj=Plot(obj,A)
            % Updating nodal coordinates after deformation with amplitude A
            obj.ncoordsplot=zeros(size(obj.ncoords));
            for i=1:obj.nodes
                obj.ncoordsplot(i,1)=obj.ncoords(i,1)+obj.u(3*i-2)*A;
                obj.ncoordsplot(i,2)=obj.ncoords(i,2)+obj.u(3*i-1)*A;
            end
            
            % Plotting
            figure; hold all; grid on; axis([min(obj.ncoords(:,1))-500 max(obj.ncoords(:,1))+1500 min(obj.ncoords(:,2))-500 max(obj.ncoords(:,2))+500]);
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2); % identify element node numbers
                
                % Plotting original structure
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2); % element node 1 - x,y original coordinates
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2); % element node 2 - x,y original coordinates
                Le=sqrt((x2-x1)^2+(y2-y1)^2); % element length
                alpha=atan2(y2-y1,x2-x1); % angle of inclination of element to the POSITIVE direction of the horizontal x axis
                plot([x1,x2],[y1,y2],'k','Linewidth',3);
                
                % Plotting the deformed structure
                x1_amp=obj.ncoordsplot(n1,1); y1_amp=obj.ncoordsplot(n1,2); % element node 1 - x,y amplified deformed coordinates
                x2_amp=obj.ncoordsplot(n2,1); y2_amp=obj.ncoordsplot(n2,2); % element node 2 - x,y amplified deformed coordinates
                % This is a little harder this time around, because the shape functions
                % are no longer linear but cubic. Recall from the lectures that:
                % the transverse deflection is v(Xi) = N1(Xi)v1 + N2(Xi)th1 + N3(Xi)v2 + N4(Xi)th2
                % the axial deflection would be u(Xi) = N5(Xi)u1 + N6(Xi)u2 where N5 and N6 are linear in Xi
                Xi=-1:0.1:+1; % The intrinsic Xi coordinate is always between -1 and +1
                
                % The following are the shape functions and 2nd derivatives evaluated at these Xi
                N1=(1/4)*(1-Xi).*(1-Xi).*(2+Xi);  ddN1ddXi=(3/2)*Xi;
                N2=(Le/8)*(1-Xi).*(1-Xi).*(1+Xi); ddN2ddXi=(1/4*(1+Xi))*Le-(1/2)*Le*(1-Xi);
                N3=(1/4)*(1+Xi).*(1+Xi).*(2-Xi);  ddN3ddXi=-(3/2)*Xi;
                N4=(Le/8)*(1+Xi).*(1+Xi).*(Xi-1); ddN4ddXi=(1/4)*Le*(-1+Xi)+(1/2*(1+Xi))*Le;
                N5=(1-Xi)/2;
                N6=(1+Xi)/2;
                
                % The following are the computed global nodal dofs
                ux1=obj.u(obj.ndofs(n1,1)); uy1=obj.u(obj.ndofs(n1,2)); thz1=obj.u(obj.ndofs(n1,3)); % element node 1 - global dofs
                ux2=obj.u(obj.ndofs(n2,1)); uy2=obj.u(obj.ndofs(n2,2)); thz2=obj.u(obj.ndofs(n2,3)); % element node 2 - global dofs
                % Transformation from global to local nodal dofs
                c=cos(alpha); s=sin(alpha);
                u1=c*ux1+s*uy1; v1=-s*ux1+c*uy1; th1=thz1; % element node 1 - local dofs
                u2=c*ux2+s*uy2; v2=-s*ux2+c*uy2; th2=thz2; % element node 1 - local dofs
                
                % The following are the v(Xi) transverse and u(Xi) axial variations within the element
                vXi=N1*v1+N2*th1+N3*v2+N4*th2;
                uXi=N5*u1+N6*u2;
                % The following is the curvature d2v/dXi2
                vppXi=ddN1ddXi*v1+ddN2ddXi*th1+ddN3ddXi*v2+ddN4ddXi*th2;
                
                % The following is a transformation of the local coordinate Xi into a
                % local element coordinate equal to 0 at its node 1 and Le_new at its node 2.
                % This is necessary because the structure is now 2D. Call this coord 's'
                s=0.5*Le*(1+Xi);
                x=x1+s*cos(alpha)-vXi*sin(alpha)*A+uXi*cos(alpha)*A;
                y=y1+s*sin(alpha)+vXi*cos(alpha)*A+uXi*sin(alpha)*A;
                
                % The tricky part - although we have just *computed* v in terms of Xi,
                % we will actually plot those same values in terms of the global x coordinate!
                for j=2:length(x)
                    % The following checks whether the element length has increased
                    % (TENSION - red) or decreased (COMPRESSION - blue)
                    tol=1e-6;
                    if u2-u1>tol
                        colourA='ro-';
                    else
                        colourA='bo-';
                    end
                    if abs(u2-u1)<tol
                        colourA='ko-';
                    end
                    
                    % The aim of computing the 2nd derivative of v(Xi) is to obtain
                    % the curvature. The bending moment is directly proportional to the
                    % curvature, so that if we know that the 2nd derivative is
                    % positive, the bending moment is also positive (SAGGING) and if
                    % the 2nd derivative is negative, the bending moment is also
                    % negative (HOGGING) - relative to how we define nodes 1 & 2
                    if vppXi(j)>=0
                        colourB='r';
                    else
                        colourB='b';
                    end
                    
                    % This ensures that the icon is plotted red for sagging and blue for hogging.
                    % Not only that, the line connecting the icons is red for axial tension and blue
                    % for axial compression. Neat!
                    plot([x(j) x(j-1)],[y(j) y(j-1)],colourA,'Linewidth',3,'Markersize',12,'MarkerFaceColor',colourB);
                end
                % This is much easier than having to change each N from Xi to x and
                % computing v in terms of x
                
                % Plotting nodes last!
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x1_amp,y1_amp,'ko','Markersize',5,'MarkerFaceColor','y');
                plot(x2_amp,y2_amp,'ko','Markersize',5,'MarkerFaceColor','y');
                
                % Display the axial load on each elements if the axial load
                % has been calculated
                if ~isempty(obj.axiload)
                    if alpha==0
                        xdisplaced=-120;
                    else
                        xdisplaced=60;
                    end
                    if abs(alpha)==pi/2
                        ydisplaced=0;
                    elseif alpha<pi/2&&alpha>0||alpha<-pi/2
                        ydisplaced=-150;
                    else
                        ydisplaced=150;
                    end
                    if abs(obj.axiload(i))<1e-6
                        text((x1_amp+x2_amp)/2+xdisplaced,(y1_amp+y2_amp)/2+ydisplaced,'0','Color','black','FontSize',12,'FontWeight','bold');
                    elseif obj.axiload(i)>0
                        text((x1_amp+x2_amp)/2+xdisplaced,(y1_amp+y2_amp)/2+ydisplaced,num2str(obj.axiload(i)),'Color','red','FontSize',12,'FontWeight','bold');
                    else
                        text((x1_amp+x2_amp)/2+xdisplaced,(y1_amp+y2_amp)/2+ydisplaced,num2str(obj.axiload(i)),'Color','blue','FontSize',12,'FontWeight','bold');
                    end
                end
            end
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);
            title(['2D frame: deformed shape amplification factor = ',num2str(A)],'FontSize',16);
        end
        
        % Display the deflection and reactions
        function obj=Display(obj)
            for dof=1:length(obj.dfree)
                disp(['The value of dof ',num2str(obj.dfree(dof)),' is ',num2str(obj.u(obj.dfree(dof)))]);
            end
            disp(' ');
            for react=1:length(obj.drestrained)
                disp(['The value of the reaction at dof ',num2str(obj.drestrained(react)),' is ',num2str(obj.f(obj.drestrained(react)))]);
            end
            disp(' ');
        end
        
        % Equilibrium checks
        function obj=EquilibriumCheck(obj)
            obj.check=zeros(1,2);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==2)))+obj.extver)<1e-6
                obj.check(2)=1;
            else
                obj.check(2)=0;
            end
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==1)))+obj.exthor)<1e-6
                obj.check(1)=1;
            else
                obj.check(1)=0;
            end
        end
        
        % Display equilibrium checks
        function obj=DisplayEquilibriumCheck(obj)
            disp('Vertical equilibrium check:');
            disp(['Total vertical reactions = ',num2str(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==2))))]);
            disp(['Total applied vertical loads = ',num2str(obj.extver)]);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==2)))+obj.extver)<1e-6
                disp('Ok.');
            end
            disp(' ');
            disp('Horizontal equilibrium check:');
            disp(['Total horizontal reactions = ',num2str(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==1))))]);
            disp(['Total applied horizontal loads = ',num2str(obj.exthor)]);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,3)==1)))+obj.exthor)<1e-6
                disp('Ok.');
            end
            disp(' ');
        end
        
        % Calculate the axial load
        function obj=AxialLoad(obj)
            obj.axiload=zeros(obj.elements,1);
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                u1=obj.u(n1*3-2); u2=obj.u(n2*3-2);
                v1=obj.u(n1*3-1); v2=obj.u(n2*3-1);
                alpha=atan2(y2-y1,x2-x1);
                c=cos(alpha); s=sin(alpha);
                Le=sqrt((x2-x1)^2+(y2-y1)^2);
                u11=u1*c+v1*s; u22=u2*c+v2*s;
                obj.axiload(i)=(obj.ea(i)/Le)*(u22-u11);
            end
        end
        
        % Buckling checks
        function obj=BucklingCheck(obj)
            % Calculate the axial load in each element (positive represents
            % tension
            obj.pcr=zeros(obj.elements,1);
            obj.mos=obj.pcr;
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                Le=sqrt((x2-x1)^2+(y2-y1)^2);
                
                % Calculate pcr
                obj.pcr(i)=obj.ei(i)*(pi^2)/((Le/2)^2);
                
                % Calculate the margin of safety
                if obj.axiload(i)<0
                    obj.mos(i)=(obj.pcr(i)/abs(obj.axiload(i)))-1;
                else
                    obj.mos(i)=NaN;
                end
            end
            
            % Calculate the minimum factor of safety
            obj.mosmin=min(obj.mos);
        end
        
        % Plot member at risk of buckling
        function obj=PlotBuckling(obj)
            figure; hold all; grid on; axis([min(obj.ncoords(:,1))-500 max(obj.ncoords(:,1))+500 min(obj.ncoords(:,2))-500 max(obj.ncoords(:,2))+500]);
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                alpha=atan2(y2-y1,x2-x1);
                if isnan(obj.mos(i))
                    plot([x1,x2],[y1,y2],'k','Linewidth',3); % black colour
                elseif obj.mos(i)>=0
                    plot([x1,x2],[y1,y2],'b','Linewidth',3); % black colour
                else
                    plot([x1,x2],[y1,y2],'r','Linewidth',3); % red colour
                end
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                
                % Display the axial load on each elements
                if alpha==0
                    xdisplaced=-100;
                else
                    xdisplaced=20;
                end
                if abs(alpha)==pi/2
                    ydisplaced=0;
                elseif alpha<pi/2&&alpha>0||alpha<-pi/2
                    ydisplaced=-90;
                else
                    ydisplaced=90;
                end
                if isnan(obj.mos(i))
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.mos(i)),'Color','black','FontSize',12,'FontWeight','bold');
                elseif obj.mos(i)>0
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.mos(i)),'Color','blue','FontSize',12,'FontWeight','bold');
                else
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.mos(i)),'Color','red','FontSize',12,'FontWeight','bold');
                end
            end
            title('Members at Risk of Buckling in Red','FontSize',20);
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);
        end
    end
end