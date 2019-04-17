%% Close all graphs and clear the memory. 

close all;
clear all;
%% Task 1
Nphotonsall = [5e4, 5e5, 5e6, 5e7];

% Calculation
% making simulated detector images corrupted by Poisson noise.

image = sum(im2double(imread('peppers.png')),3); 
norm = sum(sum(image));
image_norm = image/norm;

% loop over different photon count & noise realizations.
% modify line 22 to incorporate "imnoise"
image_noisy = zeros([size(image_norm),length(Nphotonsall)]);
for j = 1:length(Nphotonsall)
    image_noisy(:,:,j) = 1e12*imnoise(1e-12*image_norm*Nphotonsall(j),'poisson');
end

% Plot results

figure;
colormap bone;
imagesc(image_norm);

figure;
colormap bone;

subplot(2,2,1);
imagesc(image_noisy(:,:,1));
title(['image with ',num2str(Nphotonsall(1)),' photons']);
subplot(2,2,2);
imagesc(image_noisy(:,:,2));
title(['image with ',num2str(Nphotonsall(2)),' photons']);
subplot(2,2,3);
imagesc(image_noisy(:,:,3));
title(['image with ',num2str(Nphotonsall(3)),' photons']);
subplot(2,2,4);
imagesc(image_noisy(:,:,4));
title(['image with ',num2str(Nphotonsall(4)),' photons']);

%Explanation: As and how the total number of photons captured by the detector is increased, the image starts to become more clear and less noisy. Yes, the expectation holds true. 

%% Task 2
%Input of user defined parameters
% Nphotonsall = array with # of photons.
% bg = #background photons per pixel.

Nphotonsall = [1e2,2e2,5e2,1e3,2e3,5e3,1e4];
bg = 0.0;

% load images ("PSF") and ground truth positions "xtrue" and "ytrue"
% of the emitters 
load emitterPSFs PSF xtrue ytrue

% Calculation
% making arrays of pixel coordinates
Np = size(PSF,1);
Ncfg = size(PSF,3);
xy = (1:Np) - (Np+1)/2;
[Xp,Yp] = meshgrid(xy,xy);

% loop over photon count
uncerestim = zeros(length(Nphotonsall),1);
uncertheory = zeros(length(Nphotonsall),1);
for jph = 1:length(Nphotonsall)
  Nph = Nphotonsall(jph);
  
  % make images corrupted with Poisson-noise
  % dummy code now, use code of task 6.1 here
  allspots = PSF;
  
  % MLE-estimate of parameters
  % program code to find the estimates for photon count "Nest",
  % emitter position "xest" and "yest", and spot width "sigest"
  % using the given MLE-expressions.
  Nphest = zeros(Ncfg,1);
  xest = zeros(Ncfg,1);
  yest = zeros(Ncfg,1);
  sigest = zeros(Ncfg,1);
  for jcfg = 1:Ncfg
    spot = allspots(:,:,jcfg);
    image_noisy = 1e12*imnoise(1e-12*spot*Nphotonsall(jph),'poisson');
    Nphest(jcfg) = sum(sum(image_noisy));
    xest(jcfg) = 1/Nphest(jcfg)*sum(sum(image_noisy.*Xp));
    yest(jcfg) = 1/Nphest(jcfg)*sum(sum(image_noisy.*Yp));
    sigest(jcfg) = 1/(2*Nphest(jcfg))*sum(sum(image_noisy.*((Xp-xest(jcfg)).^2+(Yp-yest(jcfg)).^2)));
    sigest(jcfg) = sqrt(sigest(jcfg));
  end
  
  % uncertainty estimation
  % dummy values now, must be replaced with right expressions
  uncertheory(jph) = std(sqrt((xest-xtrue).^2+(yest-ytrue).^2));
  uncerestim(jph) = mean((sigest./sqrt(Nphest)));
end
figure;
loglog(Nphotonsall,uncertheory,'b',...
       Nphotonsall,uncerestim,'or');
axis([1e2 1e4 0.01 1.0]);
xlabel('photon count');
ylabel('localization uncertainty (pixels)');
legend('theory','simulation');
grid off;

%Explanation: There are 500 images. Each of the images have 15X15 pixels.
%The estimated localization uncertainity exhibits a better fit with the
%theoritical localization uncertainity with the increase in the number of photons. 
%% Task 3
% Input of user defined parameters
% Nphotonsall = array with # of photons.
% bg = #background photons per pixel.


Nphotonsall = [1e2,2e2,5e2,1e3,2e3,5e3,1e4];
bg = [0.0, 0.5, 3.0, 5.0, 15.0, 30.0];

% load images ("PSF") and ground truth positions "xtrue" and "ytrue"
% of the emitters 
load emitterPSFs PSF xtrue ytrue

% Calculation
% making arrays of pixel coordinates
Np = size(PSF,1);
Ncfg = size(PSF,3);
xy = (1:Np) - (Np+1)/2;
[Xp,Yp] = meshgrid(xy,xy);

% loop over photon count
uncerestim = zeros(length(Nphotonsall),1);
uncertheory = zeros(length(Nphotonsall),1);
for b = 1:length(bg)
for jph = 1:length(Nphotonsall)
  Nph = Nphotonsall(jph);
  
  % make images corrupted with Poisson-noise
  % dummy code now, use code of task 6.1 here
  allspots = PSF;
  
  % MLE-estimate of parameters
  % program code to find the estimates for photon count "Nest",
  % emitter position "xest" and "yest", and spot width "sigest"
  % using the given MLE-expressions.
  Nphest = zeros(Ncfg,1);
  xest = zeros(Ncfg,1);
  yest = zeros(Ncfg,1);
  sigest = zeros(Ncfg,1);
  for jcfg = 1:Ncfg
    spot = allspots(:,:,jcfg);
    image_noisy = 1e12*imnoise(1e-12*(spot*Nphotonsall(jph)+bg(b)),'poisson');
    Nphest(jcfg) = sum(sum(image_noisy));
    xest(jcfg) = 1/Nphest(jcfg)*sum(sum(image_noisy.*Xp));
    yest(jcfg) = 1/Nphest(jcfg)*sum(sum(image_noisy.*Yp));
    sigest(jcfg) = 1/(2*Nphest(jcfg))*sum(sum(image_noisy.*((Xp-xest(jcfg)).^2+(Yp-yest(jcfg)).^2)));
    sigest(jcfg) = sqrt(sigest(jcfg));
  end
  
  % uncertainty estimation
  % dummy values now, must be replaced with right expressions
  uncertheory(jph) = std(sqrt((xest-xtrue).^2+(yest-ytrue).^2));
  uncerestim(jph) = mean((sigest./sqrt(Nphest)));
end
figure('name',"Background = "+string(bg(b)));
loglog(Nphotonsall,uncertheory,'b',...
       Nphotonsall,uncerestim,'or');
axis([1e2 1e4 0.01 1.0]);
xlabel('photon count');
ylabel('localization uncertainty (pixels)');
legend('theory','simulation');
grid off;
end

%Explanation: With the increase in the addition of the constant background
%after 0.5 the estimated value strays away from the true value