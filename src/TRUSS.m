classdef TRUSS
    
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
        
        % rcond of Kff
        krco
        
        % Equilibrium checks
        exthor
        extver
        check
        
        % Buckling checks: Euler buckling, axial load, margin of safety
        pcr
        axiload
        mos
        mosmin
        
        % Minimum radius to prevent buckling
        rmin
    end
    
    methods
        % Creation of a variable with this class
        function obj=TRUSS(E,R,NCOORDS,ENODES,DFREE,DRESTRAINED)
            obj.e=E;
            obj.r=R;
            obj.ncoords=NCOORDS;
            obj.enodes=ENODES;
            obj.dfree=DFREE;
            obj.drestrained=DRESTRAINED;
            
            % Calculate ea and ei for each element
            obj.ea=obj.e.*(pi*(obj.r.^2));
            obj.ei=obj.e.*(pi*(obj.r.^4)/4);
            
            % Calculate the number of nodes and elements
            obj.nodes=size(obj.ncoords,1);
            obj.elements=size(obj.enodes,1);
            
            % Defining the numbering of dofs
            obj.ndofs=zeros(obj.nodes,2);
            for i=1:obj.nodes
                obj.ndofs(i,:)=[2*i-1,2*i];
            end
        end
        
        % Calculate the total volume
        function obj=Volume(obj)
            obj.volume=0;
            for i=1:obj.elements;
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
            obj.k=zeros(2*(obj.nodes));
            for i=1:obj.elements;
                n1=obj.enodes(i,1); n2=obj.enodes(i,2); % identify element node numbers
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2); % element node 1 - x,y coordinates
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2); % element node 2 - x,y coordinates
                dof11=obj.ndofs(n1,1); dof12=obj.ndofs(n1,2); % element node 1 - dofs
                dof21=obj.ndofs(n2,1); dof22=obj.ndofs(n2,2); % element node 2 - dofs
                alpha=atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE direction of the x axis
                c=cos(alpha); c2=c*c; s=sin(alpha); s2=s*s; cs=c*s; % angle parameters
                Le=sqrt((x2-x1)^2+(y2-y1)^2); % element length
                
                % Calculate ke
                ke=obj.ea(i)/Le; % element axial stiffness
                
                % Updating global stiffness matrix [K] coefficients
                % Row 1 - element dof11
                obj.k(dof11,dof11)=obj.k(dof11,dof11)+ke*c2; % Col 1 - element dof11
                obj.k(dof11,dof12)=obj.k(dof11,dof12)+ke*cs; % Col 2 - element dof12
                obj.k(dof11,dof21)=obj.k(dof11,dof21)-ke*c2; % Col 3 - element dof21
                obj.k(dof11,dof22)=obj.k(dof11,dof22)-ke*cs; % Col 4 - element dof22
                
                % Row 2 - element dof12
                obj.k(dof12,dof11)=obj.k(dof12,dof11)+ke*cs; % Col 1 - element dof11
                obj.k(dof12,dof12)=obj.k(dof12,dof12)+ke*s2; % Col 2 - element dof12
                obj.k(dof12,dof21)=obj.k(dof12,dof21)-ke*cs; % Col 3 - element dof21
                obj.k(dof12,dof22)=obj.k(dof12,dof22)-ke*s2; % Col 4 - element dof22
                
                % Row 3 - element dof21
                obj.k(dof21,dof11)=obj.k(dof21,dof11)-ke*c2; % Col 1 - element dof11
                obj.k(dof21,dof12)=obj.k(dof21,dof12)-ke*cs; % Col 2 - element dof12
                obj.k(dof21,dof21)=obj.k(dof21,dof21)+ke*c2; % Col 3 - element dof21
                obj.k(dof21,dof22)=obj.k(dof21,dof22)+ke*cs; % Col 4 - element dof22
                
                % Row 4 - element dof22
                obj.k(dof22,dof11)=obj.k(dof22,dof11)-ke*cs; % Col 1 - element dof11
                obj.k(dof22,dof12)=obj.k(dof22,dof12)-ke*s2; % Col 2 - element dof12
                obj.k(dof22,dof21)=obj.k(dof22,dof21)+ke*cs; % Col 3 - element dof21
                obj.k(dof22,dof22)=obj.k(dof22,dof22)+ke*s2; % Col 4 - element dof22
            end
            
            % Initiate U and specify restrained dof
            obj.u=zeros(2*obj.nodes,1);
            obj.u(obj.drestrained)=U;
            
            % Construct the global nodal force vector and calculate the
            % total external horizontal and vertical force to be used in
            % equilibrium checks
            obj.f=zeros(2*obj.nodes,1);
            obj.exthor=0;
            obj.extver=0;
            for i=1:size(F,1)
                obj.f(F(i,1))=F(i,2);
                if rem(F(i,1),2)==1
                    obj.exthor=obj.exthor+F(i,2);
                else
                    obj.extver=obj.extver+F(i,2);
                end
            end
            
            % Solve for the unknown nodal dofs
            obj.u(obj.dfree)=obj.k(obj.dfree,obj.dfree)\(obj.f(obj.dfree)-obj.k(obj.dfree,obj.drestrained)*obj.u(obj.drestrained));
            
            % Solve for the unknown reactions
            obj.f(obj.drestrained)=obj.k(obj.drestrained,obj.dfree)*obj.u(obj.dfree)+obj.k(obj.drestrained,obj.drestrained)*obj.u(obj.drestrained);
            
            % Updating nodal coordinates after deformation
            count=1;
            obj.ncoordsnew=zeros(size(obj.ncoords));
            for i=1:obj.nodes;
                for j=1:2;
                    obj.ncoordsnew(i,j)=obj.ncoords(i,j)+obj.u(count);
                    count=count+1;
                end
            end
        end
        
        % Find the maximum lateral sway
        function obj=MaximumLateralSway(obj)
            ux=obj.u(1:2:end);
            obj.uxmax=max(ux);
        end
        
        function obj=Condition(obj)
            % Calculate the rcond of the matrix that is inverted
            obj.krco=rcond(obj.k(obj.dfree,obj.dfree));
        end
        
        % Plot the final structure
        function obj=Plot(obj,A)
            % Updating nodal coordinates after deformation with amplitude A
            count=1;
            obj.ncoordsplot=zeros(size(obj.ncoords));
            for i=1:obj.nodes;
                for j=1:2;
                    obj.ncoordsplot(i,j)=obj.ncoords(i,j)+obj.u(count)*A;
                    count=count+1;
                end
            end
            
            % Plotting
            figure; hold all; grid on; axis([min(obj.ncoords(:,1))-500 max(obj.ncoords(:,1))+500 min(obj.ncoords(:,2))-500 max(obj.ncoords(:,2))+500]);
            for i=1:obj.elements;
                n1=obj.enodes(i,1); n2=obj.enodes(i,2); % identify element node numbers
                
                % Plotting original structure
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2); % element node 1 - x,y original coordinates
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2); % element node 2 - x,y original coordinates
                Le=sqrt((x2-x1)^2+(y2-y1)^2 ); % element length
                plot([x1,x2],[y1,y2],'k','Linewidth',3);
                
                % Check on changes in member lengths and plotting amplified deformed structure
                x1_new=obj.ncoordsnew(n1,1); y1_new=obj.ncoordsnew(n1,2); % element node 1 - x,y actual deformed coordinates
                x2_new=obj.ncoordsnew(n2,1); y2_new=obj.ncoordsnew(n2,2); % element node 2 - x,y actual deformed coordinates
                x1_amp=obj.ncoordsplot(n1,1); y1_amp=obj.ncoordsplot(n1,2); % element node 1 - x,y amplified deformed coordinates
                x2_amp=obj.ncoordsplot(n2,1); y2_amp=obj.ncoordsplot(n2,2); % element node 2 - x,y amplified deformed coordinates
                Le_new=sqrt((x2_new-x1_new)^2+(y2_new-y1_new)^2 ); % actual new element length
                if Le_new<=Le; % element length has decreased - member in compression
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                elseif Le_new>Le; % element length as increased - member in tension
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                end
                
                % Plotting nodes last!
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                
                % Display the axial load on each elements if the axial load
                % has been calculated
                alpha=atan2(y2-y1,x2-x1);
                if ~isempty(obj.axiload)
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
        end
        
        % Display the deflection and reactions
        function obj=Display(obj)
            for dof=1:length(obj.dfree);
                disp(['The value of dof ',num2str(obj.dfree(dof)),' is ',num2str(obj.u(obj.dfree(dof)))]);
            end
            disp(' ');
            for react=1:length(obj.drestrained);
                disp(['The value of the reaction at dof ',num2str(obj.drestrained(react)),' is ',num2str(obj.f(obj.drestrained(react)))]);
            end
            disp(' ');
        end
        
        % Equilibrium checks
        function obj=EquilibriumCheck(obj)
            obj.check=zeros(1,2);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==0)))+obj.extver)<1e-6
                obj.check(2)=1;
            else
                obj.check(2)=0;
            end
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==1)))+obj.exthor)<1e-6
                obj.check(1)=1;
            else
                obj.check(1)=0;
            end
        end
        
        % Display equilibrium checks
        function obj=DisplayEquilibriumCheck(obj)
            disp('Vertical equilibrium check:');
            disp(['Total vertical reactions = ',num2str(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==0))))]);
            disp(['Total applied vertical loads = ',num2str(obj.extver)]);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==0)))+obj.extver)<1e-6
                disp('Ok.');
            end
            disp(' ');
            disp('Horizontal equilibrium check:');
            disp(['Total horizontal reactions = ',num2str(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==1))))]);
            disp(['Total applied horizontal loads = ',num2str(obj.exthor)]);
            if abs(sum(obj.f(obj.drestrained(rem(obj.drestrained,2)==1)))+obj.exthor)<1e-6
                disp('Ok.');
            end
            disp(' ');
        end
        
        % Calculate the axial load in the elements
        function obj=AxialLoad(obj)
            obj.axiload=zeros(obj.elements,1);
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                u1=obj.u(n1*2-1); u2=obj.u(n2*2-1);
                v1=obj.u(n1*2); v2=obj.u(n2*2);
                alpha=atan2(y2-y1,x2-x1);
                c=cos(alpha); s=sin(alpha);
                Le=sqrt((x2-x1)^2+(y2-y1)^2);
                u11=u1*c+v1*s; u22=u2*c+v2*s;
                obj.axiload(i)=(obj.ea(i)/Le)*(u22-u11);
            end
        end
        
        % Buckling checks
        function obj=BucklingCheck(obj)
            obj.pcr=zeros(obj.elements,1);
            obj.mos=obj.pcr;
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                Le=sqrt((x2-x1)^2+(y2-y1)^2);
                
                % Calculate pcr
                obj.pcr(i)=obj.ei(i)*(pi^2)/(Le^2);
                
                % Calculate the margin of safety
                if obj.axiload(i)<0
                    obj.mos(i)=(obj.pcr(i)/abs(obj.axiload(i)))-1;
                else
                    obj.mos(i)=NaN;
                end
            end
            
            % Calculate the minimum margin of safety
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
                elseif abs(obj.mos(i))<1e-6
                    plot([x1,x2],[y1,y2],'b','Linewidth',3); % black colour
                elseif obj.mos(i)>=0
                    plot([x1,x2],[y1,y2],'b','Linewidth',3); % black colour
                else
                    plot([x1,x2],[y1,y2],'r','Linewidth',3); % red colour
                end
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                
                % Display the margin of safety on each elements
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
                elseif abs(obj.mos(i))<1e-6
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,'0','Color','blue','FontSize',12,'FontWeight','bold');
                elseif obj.mos(i)>=0
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.mos(i)),'Color','blue','FontSize',12,'FontWeight','bold');
                else
                    text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.mos(i)),'Color','red','FontSize',12,'FontWeight','bold');
                end
            end
            title('Members at Risk of Buckling in Red','FontSize',20);
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);
        end
        
        % Find the minimum radius of member at risk of buckling such that
        % mos is equal to 0
        function obj=BucklingRadius(obj)
            obj.rmin=NaN(obj.elements,1);
            for i=1:obj.elements
                if obj.mos(i)<0
                    n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                    x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                    x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                    Le=sqrt((x2-x1)^2+(y2-y1)^2);
                    obj.rmin(i)=((4*(Le^2)*abs(obj.axiload(i)))/((pi^3)*obj.e(i)))^(1/4); 
                end
            end
        end
        
        % Plot a figure showing the radius of each element
        function obj=PlotRadius(obj)
            figure; hold all; grid on; axis([min(obj.ncoords(:,1))-500 max(obj.ncoords(:,1))+500 min(obj.ncoords(:,2))-500 max(obj.ncoords(:,2))+500]);
            for i=1:obj.elements
                n1=obj.enodes(i,1); n2=obj.enodes(i,2);
                x1=obj.ncoords(n1,1); y1=obj.ncoords(n1,2);
                x2=obj.ncoords(n2,1); y2=obj.ncoords(n2,2);
                alpha=atan2(y2-y1,x2-x1);
                plot([x1,x2],[y1,y2],'k','Linewidth',3); % black colour
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                
                % Display the radius on each elements
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
                text((x1+x2)/2+xdisplaced,(y1+y2)/2+ydisplaced,num2str(obj.r(i)),'Color','black','FontSize',12,'FontWeight','bold');
            end
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);
        end
    end
end