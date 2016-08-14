if 1
    %% Reconstruction Demo
% This demo illustrates the image reconstruction capabilities of 
% OpenRecon, and is a follow-up of the Projection Demo at
% <http://www.imaging.sbes.vt.edu/software/openrecon/ OpenRecon Home> 
%
% * OpenRecon allows 2 scanning geometries --
% parallel- and fan-beam. 
% * All configuration details are specified through the use of a 
% configuration file (Config file). 
% * Furthermore there is a utility function which allows easy 
% generation of Config files in the format required by OpenRecon -- use of 
% this utility function is demonstrated. 
% * OpenRecon allows projection and reconstruction at
% non-linearly spaced angles, and this capability is demonstrated. 
% * Various reconstruction algorithms are supported: Analytic (FBP),
% Iterative (SART, SIRT, OS-SART)
% * Various image regularization schemes are supported - TV minimization,
% Soft-thresholding. The use of these are demonstrated here.
%
% <http://www.imaging.sbes.vt.edu/software/openrecon/ OpenRecon Home>

%% NOTE
% This documentation is automatically generated through execution of the
% following script in the OpenRecon source-code distribution:
%
%   demos/demo_xrec.m
%
% It is recommended that you run the same program on your local machine  
% (i.e. after installing OpenRecon),
% and compare with the results in this document. Furthermore you should
% tweak the settings in the Config file to see how your changes affect the
% results.

%% Initialization
% To use OpenRecon, we need to add the folder for OpenRecon to the MATLAB path.
% Here, we use the Shepp-Logan phantom for testing OpenRecon reconstruction
% algorithms.
clear; 
clc; 
close all;
addpath('~/Coding/recon/openrecon/');		%% replace this with your local path
addpath('~/Coding/recon/openrecon/utilities/');	%% replace this with your local path
addpath(genpath('~/Documents/Coding/utilitiesmatlab/'));	%% replace this with your local path

% displacement of the shepp logan phantom
deltaX = 5;
deltaY = 8;
% estimates of shepp logan ellipse major and minor axis
b=177;
a=233;


% FBP adjustment
rrNorm = 4;

% interior recon
rr = 0.23;      % truncation ration
sigma = 25;     % sigma of gaussian for extrapolating truncated sinogram

% plot limits
plotLim = 80;

% geometry parameters
rowN = 512; colN = rowN;
N = 360;
a0 = 0 ; %- 107;
aN = 359; %107;
d = (aN - a0)/(N-1);
theta = a0: d: aN;
pxlW = 0.9766; 
pxlH = pxlW; 
scanR = 625.61; 
detrR = 0;

% generate the displaced phantom
I1 = phantom(rowN)';
I2 = zeros(size(I1));
I2(1:end - deltaY, 1: end-deltaX) = I1(deltaY+1:end, deltaX+1:end);
I = I2;
figure; imshow(I, []);
%% generate ellipse with major and minor axis specified by scout estimates
X0=-deltaX;
Y0=-deltaY;

phi=0;
[x y] = meshgrid(-rowN/2 + 1:rowN/2,-rowN/2 + 1:rowN/2);
el=((x-X0)/a).^2+((y-Y0)/b).^2<=1;
el = el * 1;
% imagesc(imrotate(el,phi)); colormap(bone) 

I_thresh = (bwconvhull(I > 0)) * 1;

imshow(el - I_thresh, [])
diff = el - I_thresh;
nnz(diff(:))

%% Forward projection
% First we generate parallel- and fan-beam projections using the
% Shepp-Logan phantom. For information on projection generation in OpenRecon, see the
% Projection Demo at
% <http://www.imaging.sbes.vt.edu/software/openrecon/ OpenRecon Home> and
% consult the Reference Manual

% display('*******');
% display('Calculating projections in parallel beam mode');
% P1 = xprj(I, 'Rec_PARA.cfg');
% 
display('*******');
display('Calculating projections in fan beam equispaced detector mode');

fname = 'globalBP.cfg';
detN = floor(rowN*sqrt(2))+10;  % sufficient detectors to cover the image diagonal
createPrjCFGFanED(fname, I, scanR, detrR, pxlW, pxlH, ...
    a0, d, length(theta), 'detN', detN);    % the projection config file
P3 = xprj(I, fname);
createCFGFanED(fname, P3, scanR, detrR, pxlW, pxlH, ...
    'a0', a0, 'dAng', d,  ...
    'rowN', rowN, 'colN', rowN, ...
    'RecType', 'FBP');                      % the recon config file
figure; imshow(P3, [])
% P3 = xprj(I, 'Rec_FAN_ED.cfg');
% P5 = xprj(I, 'Proj_FAN_EA.cfg'); % FAN_EA scanning mode is not fully supported in OpenRecon yet


%% FBP: Fan-beam reconstruction
display('*** FBP of full-scan ****');
display('	(1) Weight for short scan');
display('	(2) Weight for fan-beam');
display('	(3) Filter projections');
display('	(4) Backproject -- parallel beam mode');
filter = 'hann'; 
P3b = P3; %weightForShortScan(P3, scanR/pxlW, (a0: d: aN)*pi/180);
P3c = weightFanProjection(P3b, scanR/pxlW, 0);
[P3d, H] = MFilterProjections(P3c, filter);
I3 = xrec(P3d, fname);
I3 = I3 / rrNorm;
%% Interior recon 1: extrapolation to fixed length
display('*** FBP of full-scan ****');
display('	(1) Weight for short scan');
display('	(2) Weight for fan-beam');
display('	(3) Filter projections');
display('	(4) Backproject -- parallel beam mode');

%% truncation of sinogram
len = size(P3b,1);
truncPtL = round(rr*len);
truncPtR = round((1-rr)*len);
N1 = truncPtL + 1;
N2 = truncPtR - 1;
P3b = P3;
P3b(1:truncPtL, :) = 0; % retain only P3(end/2-end/4: end/2+end/4, :); 
P3b(truncPtR:end, :) = 0; 
figure; imshow(P3b, []);

%% extrapolation of sinogram by 'edgegauss' method (from "extendSino.m" in
% utilitiesmatlab repo
% sinogram extension by pulling out the border values
ESino = P3b;
for j = 1: size(ESino, 2) %for all angles
    endval1 = ESino(N1, j); 
    if endval1 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        for i = 1: N1 - 1
            ESino(i, j) = endval1;
        end
    end

    endval2 = ESino(N2, j); 
    if endval2 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        for i = N2 + 1: size(ESino, 1)
            ESino(i, j) = endval2;
        end
    end
end
imshow(ESino, []);

%
% apply edge Gauss to sinogram
PR2g = ESino;
for j = 1: size(PR2g, 2)
    arr = PR2g(:, j);
    arr = edgeGauss1D_v2(arr, N1, sigma);
    PR2g(:, j) = arr;
end
P3c = PR2g;
figure; imshow(PR2g, [])

%

P3d = weightFanProjection(P3c, scanR/pxlW, 0);
[P3e, H] = MFilterProjections(P3d, filter);
I4 = xrec(P3e, fname);
I4 = I4 / rrNorm;


%% generate projections for estimated ellipse
pHull = xprj(el, fname);
%% extrapolation of sinogram by 'edgegauss' to edge of ellipse projection
ESino = P3b;
for j = 1: size(ESino, 2) %for all angles
    endval1 = ESino(truncPtL + 1, j); 
    if endval1 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        for i = 1: N1 - 1
            ESino(i, j) = endval1;
        end
    end

    endval2 = ESino(truncPtR - 1, j); 
    if endval2 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        for i = N2 + 1: size(ESino, 1)
            ESino(i, j) = endval2;
        end
    end
end
figure; imshow(ESino, []);

%
% apply edge Gauss to sinogram
PR2g = ESino;
for j = 1: size(PR2g, 2)
    arr = PR2g(:, j);
    xHull = pHull(:, j);
    
    endval1 = ESino(N1, j); 
    if endval1 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        N1_hull = find(xHull, 1, 'first' );
        if (N1 > N1_hull)
            sigmaL = (N1-N1_hull) / 1;
            arr = edgeGauss1D_left(arr, N1, sigmaL);
        else
            arr(1: truncPtL) = 0;
        end
    end
    
    endval2 = ESino(N2, j); 
    if endval2 ~= 0 %if non-zero, pull down to zero. Otherwise do not need to do anything
        N2_hull = find(xHull, 1, 'last' );
        if (N2_hull > N2)
            sigmaR = (N2_hull-N2) / 1;
            arr = edgeGauss1D_right(arr, N1, sigmaR);
        else
            arr(truncPtR: end) = 0;
        end
    end
    PR2g(:, j) = arr;
end
P3c = PR2g;
figure; imshow(PR2g, [])

P3d = weightFanProjection(P3c, scanR/pxlW, 0);
[P3e, H] = MFilterProjections(P3d, filter);
I5 = xrec(P3e, fname);
I5 = I5 / rrNorm;




%% Figure displays
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)*0.3]),
subplot(2,2,1), imshow(I, []); title('Shepp-Logan phantom');
% subplot(1,3,2), imshow(I1, []); title('Parallel-beam recon');
subplot(2,2,2), imshow(I3, []); title('Fan-beam recon');
subplot(2,2,3), imshow(I4, []); title('Edge-gauss extrapolation');
subplot(2,2,4), imshow(I5, []); title('Scout based');

%%
figure, 
imshow([[I, I3]; [I4, I5]], []);
%%
nrow = 167;
ncol = [102 403];
figure;
subplot(2,1,1)
% plot( I(end/2, end/2-plotLim: end/2+plotLim))
% hold on
% plot(I3(end/2, end/2-plotLim: end/2+plotLim), 'r')
% plot(I4(end/2, end/2-plotLim: end/2+plotLim), 'k--')
% plot(I5(end/2, end/2-plotLim: end/2+plotLim), 'b--')
plot( I(nrow, ncol(1): ncol(2) ))
hold on
plot(I3(nrow, ncol(1): ncol(2)), 'r')
plot(I4(nrow, ncol(1): ncol(2)), 'k--')
plot(I5(nrow, ncol(1): ncol(2)), 'b--')
hold off;
legend('Phantom', 'FBP', 'Edge-gauss', 'Scout IT');

subplot(2,1,2)
plot(I(end/2-plotLim: end/2+plotLim, end/2))
hold on
plot(I3(end/2-plotLim: end/2+plotLim, end/2), 'r')
plot(I4(end/2-plotLim: end/2+plotLim, end/2), 'k--')
plot(I5(end/2-plotLim: end/2+plotLim, end/2), 'b--')
hold off;
legend('Phantom', 'FBP', 'Edge-gauss', 'Scout IT');

end

%% modified scout configuration
dd = 200;
scanR2 = (scanR + dd);
detrR2 = 0;
pxlW2 = (scanR + dd)/scanR * pxlW;
pxlH2 = pxlW2;
% rrAP = (scanR + dd)/scanR;
fname = 'globalBP2.cfg';
% detN = floor(rowN*sqrt(2))+10;  % sufficient detectors to cover the image diagonal
I2 = imresize(I, pxlW/pxlW2);
createPrjCFGFanED(fname, I2, scanR2, detrR2, pxlW2, pxlH2, ...
    a0, d, length(theta), 'detN', detN);    % the projection config file
P3_AP = xprj(I2, fname);
figure; 
hold on
plot(P3(:,1));
plot(P3_AP(:,1), '.')
legend('AP orig', 'AP mod');
%% visualize the original AP scout (if there were no truncation) and in modified configuration
D = 1097.61;        % Source to detector distance
xIso =(-len/2+0.5:len/2-0.5) * pxlW;
xAP = (-len/2+0.5:len/2-0.5) * pxlW2;
xDet = (-len/2+0.5:len/2-0.5) * D/scanR * pxlW;
% P3x = interp1(xIso, P3, xDet);
% P3_APx = interp1(xAP, P3_AP, xDet);
%%
figure;
plot(xDet, P3(:,1))
hold on
plot(xDet, P3_AP(:,1))
plot(xDet, P3(:, 90))
hold off;
legend('AP orig', 'AP mod', 'ML');

%% Get results from graph and get the set of equations
% Geometry
r0 = scanR;         % Source to isocenter distance
% D = 1097.61;        % Source to detector distance
% pxlW = 0.625;

% Experiments
% dd = 625  ;            % displacement of the detector in AP scout
% nDet = size(P3,1);           % number of detectors in simulation

% Recorded observations
% iAP = [120.5, 603.5];
% iML = [188, 565.5];
% % Transformation to mm values
% p12 = (iAP - nDet/2) * pxlW * D/r0;
% p34 = (iML - nDet/2) * pxlW * D/r0;
% p1 = p12(1);    p2 = p12(2);
% p3 = p34(1);    p4 = p34(2);
% 
p1 = xDet(find(P3_AP(:,1),1,'first'));
p2 = xDet(find(P3_AP(:,1),1,'last'));
p3 = xDet(find(P3(:,90),1,'first'));
p4 = xDet(find(P3(:,90),1,'last'));
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

%% Plug results into Mathematica and get the solutions

% BB = 9.97216e-5;            % Desired A of Mathematica: 1.93147e-5
% AA = 3.319e-5;              % Use B  of Mathematica
% x0Det = 5.11929;            % Keep same sign as Mathematica (but for y)
% y0Det = -3.79241;           % Desired: 8.1917
% RyDet = 1/sqrt(AA)/pxlW;
% RxDet = 1/sqrt(BB)/pxlW;
% X0Det = x0Det * pxlW;
% Y0Det = y0Det * pxlW;
RyDet = 173;
RxDet = 233;
X0Det = -8.15581;
Y0Det = -4.8863;
el=((x-X0Det)/RxDet).^2+((y-Y0Det)/RyDet).^2<=1;
el = el * 1;
I_thresh = (bwconvhull(I > 0)) * 1;
figure; imshow(el - I_thresh, [])