function f=TrussSway(R,E,NCOORDS,ENODES,DFREE,DRESTRAINED,U,F)

% Ho Zhu Yi
% 01082991
% 9 March 2017

% All dimensions in N and mm.

% Solve for deflection
truss1=TRUSS(E,R,NCOORDS,ENODES,DFREE,DRESTRAINED);
truss1=truss1.Solve(U,F);
truss1=truss1.MaximumLateralSway();

% Determine the lateral sway at the point of application of the top-most
% load
sway=truss1.uxmax;

% Minus 25 from sway to allow the use of fzero
f=sway-25;

end