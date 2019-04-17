%% Imaging systems 2017-2018 (AP3121 D), MATLAB assignment-tomography
% made by J. Kalkman, send your comments to j.kalkman@tudelft.nl
% assemble plots in single figure and put labels/units on all your axes
% upload your answer to Brightspace assignments
clc
clear all; close all;
FigHandle1 = figure(1); set(FigHandle1, 'Position', [700, 50, 1200, 900]); set(FigHandle1,'numbertitle','off','name','Filtered Back Projection', 'color','w'); 

% Define parameters
Nviews=32;                                                        % # of angular views
Nxy=32;                                                           % number of points in X and Y
D=1.0;                                                            % object dimension space [mm]
k=64;                                                             % 2k - 1 is the number of detector elements
tau=0.02;                                                         % detector pitch
t=tau*(-(k-1):(k-1))';                                            % coordinates of detector array, make sure detector array is larger than 2*D
detectorwidth=max(t)-min(t);                                      
%% 
% Create two 2D grids with X and Y coordinates of all pixels in the sample            
deltat=t(2)-t(1);
deltax=D/Nxy;
x=-D/2:deltax:D/2-deltax;
y=x;
[X,Y] = meshgrid(x,x);
f = -1/2/deltat:1/detectorwidth:1/2/deltat+(1/detectorwidth);
[Fx,Fy] = meshgrid(f,f);
corvec = exp(1i*pi*(f)*detectorwidth);
corvec2=exp(1i*pi*(Fx+Fy)*detectorwidth);
%%
% Preallocation of reconstruction arrays
BackProj_P=zeros(size(X));                                      % Backprojection
BackProj_Q=zeros(size(X));                                      % Filtered backprojection 
obj=(X.^2+Y.^2)<(D^2);
X_vec=X(obj);
Y_vec=Y(obj);
%%
object=zeros(Nxy, Nxy);
R=0.125*D;
cirob=(X.^2+Y.^2)<(R^2);
circent=cirob(12:21,12:21);
object(8:17,8:17)=circent; % shift the obj away from centre by 4 units
%%
% Plot the object
figure(1); subplot(331); contourf(x,y,object, 'edgecolor', 'none'); xlabel('x-coordinate (mm)'); ylabel('y-coordinate (mm)'); 
title('Phantom'); colormap(gray); grid on; colorbar; set(gca, 'Clim', [-0.1 1.1]); 
x0=x(12)+0.5*(D/(Nxy-1));
y0=x0;
rho=1;
theta=0;
%% Backprojection without filter
theta_vec=0:(2*pi)/(Nviews-1):2*pi;
Ptheta_t_array=zeros(2*k-1, Nviews);   % Array to be filled with projections (aka the sinogram) 
obj_pr=zeros(127,127);
object_pr=zeros(127,127);
figure

for i=1:32
    xp=x*cos(theta_vec(i))+y*sin(theta_vec(i));
    yp=-x*sin(theta_vec(i))+y*cos(theta_vec(i));
    x0p=x0*cos(theta_vec(i))+y0*sin(theta_vec(i));
    aa=(x0p-t).^2;
    pr=2*rho*sqrt((R^2)-((x0p-t).^2));
    c=aa<R^2;
    pr=c.*pr;
    Ptheta_t_array(:,i)=pr;
    Ptheta_t_array1(:,i) = corvec.*fftshift(fft(pr'))*deltat^2;
    obj_pr(64,:)=Ptheta_t_array1(:,i);
    angle=rad2deg(theta_vec(i));
    rot=imrotate(obj_pr,angle,'crop');
    object_pr=object_pr+rot;
    subplot(8,4,i)
    plot(f,real(Ptheta_t_array(:,i)))
end
backpr=corvec2.*ifftshift(ifft2(object_pr))/deltat^2;
figure
colormap(gray);
imagesc(real(backpr))
 %% Backprojection with ramp filter
Qtheta_t_array=zeros(2*k-1, Nviews); 
obj_pr=zeros(127,127);
object_pr=zeros(127,127);
figure
for i=1:32
    xp=x*cos(theta_vec(i))+y*sin(theta_vec(i));
    yp=-x*sin(theta_vec(i))+y*cos(theta_vec(i));
    x0p=x0*cos(theta_vec(i))+y0*sin(theta_vec(i));
    aa=(x0p-t).^2;
    pr=2*rho*sqrt((R^2)-((x0p-t).^2));
    c=aa<R^2;
    pr=c.*pr;
    Qtheta_t_array(:,i) = corvec.*fftshift(fft(pr'))*deltat^2;
    obj_pr(64,:)=Qtheta_t_array(:,i);
    angle=rad2deg(theta_vec(i));
    rot=imrotate(obj_pr,angle,'crop');
    object_pr=object_pr+rot;
    subplot(8,4,i)
    plot(f,real(Qtheta_t_array(:,i)))
end

H=abs(f); % Ramp filter
filter=object_pr.*H;
figure
colormap(gray);
imagesc(real(filter))

backprf=corvec2.*ifftshift(ifft2(filter))/deltat^2;
figure
colormap(gray);
imagesc(real(backprf))

% more streaks can be observed with reduced number or projections
% My observation: Number of streaks observed is higher after filtering
% because all the artifacts are covered in the background in the image that
% is not filtered. 