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
addpath('~/Documents/Coding/openrecon/');		%% replace this with your local path
addpath('~/Documents/Coding/openrecon/utilities/');	%% replace this with your local path
addpath(genpath('~/Documents/Coding/utilitiesmatlab/'));	%% replace this with your local path

N = 256;
I = phantom(N)';

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
P3 = xprj(I, 'Rec_FAN_ED.cfg');
% P5 = xprj(I, 'Proj_FAN_EA.cfg'); % FAN_EA scanning mode is not fully supported in OpenRecon yet


%% FBP: Fan-beam reconstruction
display('*** FBP of full-scan ****');
display('	(1) Weight for short scan');
display('	(2) Weight for fan-beam');
display('	(3) Filter projections');
display('	(4) Backproject -- parallel beam mode');
filter = 'hann'; 
scanR = 571.125; pxlW = 1;
a0 = -90; d = 1; numAng = 360; aN = a0 + (numAng-1)*d;
P3b = P3; %weightForShortScan(P3, scanR/pxlW, (a0: d: aN)*pi/180);
P3c = weightFanProjection(P3b, scanR/pxlW, 0);
[P3d, H] = MFilterProjections(P3c, filter);
I3 = xrec(P3d, 'Rec_FAN_ED.cfg');

%% Interior recon 1: extrapolation to fixed length
display('*** FBP of full-scan ****');
display('	(1) Weight for short scan');
display('	(2) Weight for fan-beam');
display('	(3) Filter projections');
display('	(4) Backproject -- parallel beam mode');
filter = 'hann'; 
scanR = 571.125; pxlW = 1;
a0 = -90; d = 1; numAng = 360; aN = a0 + (numAng-1)*d;

%% truncation of sinogram
len = size(P3b,1);
rr = 0.28; 
sigma = 25;
truncPtL = round(rr*len);
truncPtR = round((1-rr)*len);
N1 = truncPtL + 1;
N2 = truncPtR - 1;
P3b = P3;
P3b(1:truncPtL, :) = 0; % retain only P3(end/2-end/4: end/2+end/4, :); 
P3b(truncPtR:end, :) = 0; 
imshow(P3b, []);

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
    arr = edgeGauss1D(arr, N1, sigma);
    PR2g(:, j) = arr;
end
P3c = PR2g;
figure; imshow(PR2g, [])

%

P3d = weightFanProjection(P3c, scanR/pxlW, 0);
[P3e, H] = MFilterProjections(P3d, filter);
I4 = xrec(P3e, 'Rec_FAN_ED.cfg');


%% generate ellipse with major and minor axis specified by scout estimates
X0=0;
Y0=0;

a=233/2;
b=177/2;
phi=0;
[x y] = meshgrid(-N/2 + 1:N/2,-N/2 + 1:N/2);
el=((x-X0)/a).^2+((y-Y0)/b).^2<=1;
el = el * 1;
% imagesc(imrotate(el,phi)); colormap(bone) 

I_thresh = (bwconvhull(I > 0)) * 1;

imshow(el - I_thresh, [])
diff = el - I_thresh;
nnz(diff(:))

%% generate projections for estimated ellipse
pHull = xprj(el, 'Rec_FAN_ED.cfg');
end
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
I5 = xrec(P3e, 'Rec_FAN_ED.cfg');




%% Figure displays
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)*0.3]),
subplot(2,3,1), imshow(I, []); title('Shepp-Logan phantom');
% subplot(1,3,2), imshow(I1, []); title('Parallel-beam recon');
subplot(2,3,2), imshow(I3, []); title('Fan-beam recon');
subplot(2,3,3), imshow(I4, []); title('Fan-beam recon');
subplot(2,3,4), imshow(I5, []); title('Fan-beam recon');

figure;
subplot(2,1,1)
plot(I(end/2, end/2-80: end/2+80))
hold on
plot(I3(end/2, end/2-80: end/2+80)/4, 'r')
plot(I4(end/2, end/2-80: end/2+80)/4, 'k--')
plot(I5(end/2, end/2-80: end/2+80)/4, 'b--')
hold off;

subplot(2,1,2)
plot(I(end/2-80: end/2+80, end/2))
hold on
plot(I3(end/2-80: end/2+80, end/2)/4, 'r')
plot(I4(end/2-80: end/2+80, end/2)/4, 'k--')
plot(I5(end/2-80: end/2+80, end/2)/4, 'b--')
hold off;
