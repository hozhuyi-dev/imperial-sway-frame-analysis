function f=FrameSway(r,E,NCOORDS,ENODES,DFREE,DRESTRAINED,U,F)

% Ho Zhu Yi
% 01082991
% 9 March 2017

% All dimensions in N and mm.

% Solve for deflection
frame1=FRAME(E,r,NCOORDS,ENODES,DFREE,DRESTRAINED);
frame1=frame1.Solve(U,F);
frame1=frame1.MaximumLateralSway();

% Determine the lateral sway at the point of application of the top-most
% load
sway=frame1.uxmax;

% Minus 25 from sway to allow the use of fzero
f=sway-25;

end