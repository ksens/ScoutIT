clc;
% Geometry
r0 = scanR;         % Source to isocenter distance
% D = 1097.61;        % Source to detector distance
% pxlW = 0.625;

% Experiments
dd = 0  ;            % displacement of the detector in AP scout
nDet = size(P3,1);           % number of detectors in simulation

% Recorded observations
% iAP = [120.5, 603.5];
% iML = [188, 565.5];
% % Transformation to mm values
% p12 = (iAP - nDet/2) * pxlW * D/r0;
% p34 = (iML - nDet/2) * pxlW * D/r0;
% p1 = p12(1);    p2 = p12(2);
% p3 = p34(1);    p4 = p34(2);
% 
p1 = -241.7;
p2 = 228.8;
p3 = -174.3;
p4 = 193.9;
% Calculations for system of nonlinear equations
S1 = -D*(1/p1 + 1/p2);
Delta1 = D^2*(1/p1-1/p2)^2;
S2 = 1/D*(p3+p4);
Delta2 = 1/D^2*(p3-p4)^2;

% % Double check
% m1 = -D/p1;
% m2 = -D/p2;
% m3 = p3/D;
% m4 = p4/D;
% m1 + m2
% (m1-m2)^2
% m3 + m4
% (m3-m4)^2
% S*(B*x^2-1)+2*B*x*(r+dd-y)==0,
eqn1 = sprintf('%f*(B*x^2-1)+2*B*x*(%f + %f -y)==0,', S1, r0, dd);
% A*D*(B*x^2-1)^2-4*B*(A*(r+dd-y)^2 + B*x^2 - 1) == 0, 
eqn2 = sprintf('A * %f *(B*x^2-1)^2-4*B*(A*(%f + %f -y)^2 + B*x^2 - 1) == 0,' , Delta1 , r0, dd);
% T*(B*(x+r)^2-1)-2*B*y*(x+r)==0, 
eqn3 = sprintf('%f *(B*(x + %f)^2-1)-2*B*y*(x + %f)==0, ', S2, r0, r0);
% A*F*(B*(x+r)^2-1)^2-4*B*(A*y^2 + B*(x+r)^2 - 1) == 0
eqn4 = sprintf('A * %f *(B*(x + %f)^2-1)^2-4*B*(A*y^2 + B*(x + %f)^2 - 1) == 0', Delta2, r0, r0);

disp(eqn1);
disp(eqn2);
disp(eqn3);
disp(eqn4);

% ea = sprintf('y - %f * x - %f == 0, ', m1, r0);
% eb = sprintf('y - %f * x - %f == 0, ', m2, r0);
% ec = 
