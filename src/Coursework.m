% Coursework - Modelling of a sway structure

% Ho Zhu Yi
% 01082991
% 16 March 2017

% All dimensions in N and mm.

close all;
clear;
clc;

%% Initialization

% Define the number of bays and storeys
BAYS=3;
STOREYS=5;

% Define E and R
E=2e5;
R=25;

% Determine the coordinates of the nodes
NCOORDS=zeros((BAYS+1)*(STOREYS+1),2);
j=0;
for i=1:((BAYS+1)*(STOREYS+1))
    NCOORDS(i,:)=[(i-j*(BAYS+1)-1)*1000,j*1000];
    if rem(i,(BAYS+1))==0
        j=j+1;
    end
end

% Define base structure configurations
FLOOR=zeros((2*BAYS)-1,2);
for i=1:(BAYS+1)
    FLOOR(i,:)=[i,i+(BAYS+1)];
end
for i=(BAYS+2):(2*(BAYS+1)-1)
    FLOOR(i,:)=[i,i+1];
end
BASE=zeros((2*BAYS+1)*STOREYS,2);
for i=1:STOREYS
    BASE(((i-1)*(2*BAYS+1)+1):(i*(2*BAYS+1)),:)=FLOOR+((i-1)*(BAYS+1));
end
clear('FLOOR');

% Determine DFREE, DRESTRAINED
TDFREE=(2*(BAYS+1))+1:2*(BAYS+1)*(STOREYS+1);
TDRESTRAINED=1:(2*(BAYS+1));

% Define the restrained u
TU=zeros((2*(BAYS+1)),1);

% Define the force on the structure f (1st column if the dof, 2nd column is
% the magnitude of the force)
TF=zeros(STOREYS,2);
TF(:,2)=100000;
for i=1:STOREYS
    TF(i,1)=i*(2*(BAYS+1))+1;
end

% Determine DFREE, DRESTRAINED
FDFREE=(3*(BAYS+1))+1:3*(BAYS+1)*(STOREYS+1);
FDRESTRAINED=1:(3*(BAYS+1));

% Define the restrained u
FU=zeros((3*(BAYS+1)),1);

% Define the force on the structure f (1st column is the dof, 2nd column is
% the magnitude of the force)
FF=zeros(STOREYS,2);
FF(:,2)=100000;
for i=1:STOREYS
    FF(i,1)=i*(3*(BAYS+1))+1;
end

% Define the different possible configurations for each storey
config=zeros(2*BAYS,2);
for i=1:BAYS
    config(2*i-1,:)=[i,i+BAYS+2];
    config(2*i,:)=[i+1,i+BAYS+1];
end

%% Question 1

% Initialize the vectors for determinant and rcond values
rco1=zeros(STOREYS+1,1);

% Starting with no cross-members, add a cross-member to each storey with
% each for loop. Determine the determinant and rcond of KFF in each loop.
% Plot in a figure with subplots
ENODES=BASE;
for i=0:STOREYS
    truss1=TRUSS(E*ones((2*BAYS+1)*STOREYS+i,1),R*ones((2*BAYS+1)*STOREYS+i,1),NCOORDS,ENODES,TDFREE,TDRESTRAINED);
    truss1=truss1.Solve(TU,TF);
    truss1=truss1.Condition();
    truss1=truss1.Plot(1);
    rco1(i+1)=truss1.krco;
    ENODES=[ENODES;config(1,:)+i*(BAYS+1)];
end

% To show the importance of ONE cross-member per storey, and that the
% stability does not depend solely on the redundancy, add 5 cross-members
% to the structure but not one on each storey.
ENODES=[BASE;zeros(STOREYS,2)];
for j=1:STOREYS-1
    ENODES(end-(STOREYS-j),:)=config(1,:)+(j-1)*(BAYS+1);
end
ENODES(end,:)=config(end,:)+(j-1)*(BAYS+1);
truss1=TRUSS(E*ones((2*BAYS+1)*STOREYS+STOREYS,1),R*ones((2*BAYS+1)*STOREYS+STOREYS,1),NCOORDS,ENODES,TDFREE,TDRESTRAINED);
truss1=truss1.Solve(TU,TF);
truss1=truss1.Condition();
truss1=truss1.Plot(1);

% Display the determinant and rcond
disp('QUESTION 1 ======================================================');
disp(' ');
fprintf('rcond of K is %g. \n\n',truss1.krco);

%% Question 2 and 3
disp('QUESTION 2 AND 3 ================================================');
disp(' ');

% Initiate a matrix to store the equilibrium checks
equi=zeros((2*BAYS)^STOREYS,2);

% Initiate a vector to store lateral sway for each configuration (there is
% a total of 6^5=7776 poosible arrangements)
sway=zeros((2*BAYS)^STOREYS,1);

for i=1:(2*BAYS)^STOREYS
    % Current arrangement
    arra=ones(1,STOREYS);
    j=str2double(dec2base(i-1,2*BAYS));
    for k=1:STOREYS
        l=floor(j/(10^(STOREYS-k)));
        arra(k)=arra(k)+l;
        j=j-l*(10^(STOREYS-k));
    end
    
    % Define the structure configuration
    ENODES=[BASE;zeros(STOREYS,2)];
    for j=1:STOREYS
        ENODES(end-(STOREYS-j),:)=config(arra(j),:)+(j-1)*(BAYS+1);
    end
    
    % Solve the structure
    truss1=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),R*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES,TDFREE,TDRESTRAINED);
    truss1=truss1.Solve(TU,TF);
    truss1=truss1.MaximumLateralSway();
    truss1=truss1.EquilibriumCheck();
    
    % Determine the lateral sway in the highest level
    sway(i)=truss1.uxmax;
    
    % Equilibrium check
    equi(i,:)=truss1.check;
end

% Equilibrium check
disp('Equilibrium check: ');
if isempty(find(equi==0, 1))
    disp('Ok.');
else
    disp('Not ok.');
end
disp(' ');

% Determine the minimum and maximum lateral sway
[swaymin,i1]=min(sway);
[swaymax,i2]=max(sway);

% Display the minimum and maximum lateral sway
disp('Minimum lateral sway: ');
fprintf('%gmm\n\n',swaymin);
disp('Maximum lateral sway: ');
fprintf('%gmm\n\n',swaymax);

% Determine the arrangement from i1 and i2
% Initiate a vector arra
arra1=ones(1,STOREYS);
arra2=arra1;
j=str2double(dec2base(i1-1,2*BAYS));
for k=1:STOREYS
    l=floor(j/(10^(STOREYS-k)));
    arra1(k)=arra1(k)+l;
    j=j-l*(10^(STOREYS-k));
end
j=str2double(dec2base(i2-1,2*BAYS));
for k=1:STOREYS
    l=floor(j/(10^(STOREYS-k)));
    arra2(k)=arra2(k)+l;
    j=j-l*(10^(STOREYS-k));
end

% Display the optimal arrangement and the least optimal arrangement
disp('Optimal arrangement: ');
disp(arra1);
disp('Least optimal arrangement: ');
disp(arra2);

% Plot both the configurations
ENODES1=[BASE;zeros(STOREYS,2)];
for j=1:STOREYS
    ENODES1(end-(STOREYS-j),:)=config(arra1(j),:)+(j-1)*(BAYS+1);
end
truss1=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),R*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES1,TDFREE,TDRESTRAINED);
truss1=truss1.Solve(TU,TF);
truss1=truss1.AxialLoad();
truss1=truss1.Plot(3);
truss1=truss1.Volume();
ENODES2=[BASE;zeros(STOREYS,2)];
for j=1:STOREYS
    ENODES2(end-(STOREYS-j),:)=config(arra2(j),:)+(j-1)*(BAYS+1);
end
truss2=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),R*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES2,TDFREE,TDRESTRAINED);
truss2=truss2.Solve(TU,TF);
truss2=truss2.AxialLoad();
truss2=truss2.Plot(3);

% Display the volumes
fprintf('The volume in a truss structure is %g. \n\n',truss1.volume);

%% Question 4
disp('QUESTION 4 ======================================================');
disp(' ');

% Create the FRAME
frame1=FRAME(E*ones((2*BAYS+1)*STOREYS,1),R*ones((2*BAYS+1)*STOREYS,1),NCOORDS,BASE,FDFREE,FDRESTRAINED);
frame1=frame1.Solve(FU,FF);
frame1=frame1.MaximumLateralSway();
frame1=frame1.AxialLoad();
frame1=frame1.Plot(1);
frame1=frame1.DisplayEquilibriumCheck();
frame1=frame1.Volume();

% Display the lateral sway
fprintf('The lateral sway is %gmm\n',frame1.uxmax);
fprintf('The volume in a frame structure is %g. \n\n',frame1.volume);

%% Question 5
disp('QUESTION 5 ======================================================');
disp(' ');

% Determine r for optimal arrangement truss
f1=@(r)(TrussSway(r*ones((2*BAYS+2)*STOREYS,1),E*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES1,TDFREE,TDRESTRAINED,TU,TF));
r1=fzero(f1,25);

% Determine r for least optimal arrangement truss
f2=@(r)(TrussSway(r*ones((2*BAYS+2)*STOREYS,1),E*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES2,TDFREE,TDRESTRAINED,TU,TF));
r2=fzero(f2,25);

% Frame
f3=@(r)(FrameSway(r*ones((2*BAYS+1)*STOREYS,1),E*ones((2*BAYS+1)*STOREYS,1),NCOORDS,BASE,FDFREE,FDRESTRAINED,FU,FF));
r3=fzero(f3,25);

% Display
disp('The minimum cross-sectional radius: ');
fprintf('For optimal arrangement truss is %gmm. \n',r1);
fprintf('For least-optimal arrangement truss is %gmm. \n',r2);
fprintf('For frame is %gmm. \n',r3);
disp(' ');

% Determine volume for optimal arrangement truss
truss1=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),r1*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES1,TDFREE,TDRESTRAINED);
truss1=truss1.Volume();
v1=truss1.volume;

% Determine volume for least-optimal arrangement truss
truss2=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),r2*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES2,TDFREE,TDRESTRAINED);
truss2=truss2.Volume();
v2=truss2.volume;

% Determine volume for frame
frame1=FRAME(E*ones((2*BAYS+1)*STOREYS,1),r3*ones((2*BAYS+1)*STOREYS,1),NCOORDS,BASE,FDFREE,FDRESTRAINED);
frame1=frame1.Volume();
v3=frame1.volume;

% Display the volumes
disp('The volume: ');
fprintf('For optimal arrangement truss is %gm^3. \n',v1/(10^9));
fprintf('For least-optimal arrangement truss is %gm^3. \n',v2/(10^9));
fprintf('For frame is %gm^3. \n',v3/(10^9));
disp(' ');

%% Question 6
disp('QUESTION 6 ======================================================');
disp(' ');

% Optimal arrangement truss
truss1=truss1.Solve(TU,TF);
truss1=truss1.AxialLoad();
truss1=truss1.BucklingCheck();
truss1=truss1.PlotBuckling();

% Least-optimal arrangement truss
truss2=truss2.Solve(TU,TF);
truss2=truss2.AxialLoad();
truss2=truss2.BucklingCheck();
truss2=truss2.PlotBuckling();

% Frame
frame1=frame1.Solve(FU,FF);
frame1=frame1.AxialLoad();
frame1=frame1.BucklingCheck();
frame1=frame1.PlotBuckling();

% Display the minimum factor of safety for each structure
disp('The margin of safety: ');
fprintf('For optimal arrangement truss is %g. \n',truss1.mosmin);
fprintf('For least-optimal arrangement truss is %g. \n',truss2.mosmin);
fprintf('For frame is %g. \n',frame1.mosmin);
disp(' ');

%% Further optimisation
disp('FURTHER OPTIMISATION ============================================');
disp(' ');

% Find the minimum radius such that mos=0
truss1=truss1.BucklingRadius(); 
R4=truss1.rmin; 

% Find the radius of other members
f4=@(r)(TrussOptimisation(r,R4,E*ones((2*BAYS+2)*STOREYS,1),NCOORDS,ENODES1,TDFREE,TDRESTRAINED,TU,TF));
r4=fzero(f4,25);
R4(isnan(R4))=r4; 

% Optimised optimal arrangement truss
truss3=TRUSS(E*ones((2*BAYS+2)*STOREYS,1),R4,NCOORDS,ENODES1,TDFREE,TDRESTRAINED);
truss3=truss3.Volume();
v4=truss3.volume;
truss3=truss3.Solve(TU,TF);
truss3=truss3.AxialLoad();
truss3=truss3.BucklingCheck();
truss3=truss3.PlotBuckling();
truss3=truss3.PlotRadius(); 

% Display the volume
fprintf('The volume for the optimised optimal arrangement truss is %gm^3. \n',v4/(10^9));