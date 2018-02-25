% ANISODIFF - Anisotropic diffusion.
%
% 
%  diff = anisodiff(im, niter, kappa, lambda, option)
%
% 
%         im     - input image
%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2
%                  3 MAD - Median Absolute Deviation (Black et el.)
%
% Return
%         diff   - diffused image.
%
% kappa controls conduction as a function of gradient.  If kappa is low
% then mall intensity gradients are able to block conduction and hence diffusion
% across step edges.  A large value reduces the influence of intensity
% gradients on conduction.
%
% lambda controls speed of diffusion (you usually want it at a maximum of
% 0.25)
%
% Diffusion equation 1 preserve high contrast edges over low contrast ones.
% Diffusion equation 2 favours wide regions over smaller ones.
% Edit : Diffusion coeffecients can be computed using Median Absolute Deviation
% value [references : 
% Michael J. Black, Member, IEEE, Guillermo Sapiro, Member, IEEE,
% David H. Marimont, Member, IEEE, and David Heeger
% Robust Anisotropic Diffusion
% IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 7, NO. 3, MARCH 1998]

% Reference: 
% P. Perona and J. Malik. 
% Scale-space and edge detection using anisotropic diffusion.
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% 12(7):629-639, July 1990.


function diff = anisodiff(im, niter, kappa, lambda, option,im_good)

if ndims(im)==3
  im = rgb2gray(im);
  im_good = rgb2gray(im_good);
end

im = double(im);
im_good = double(im_good);
[rows,cols] = size(im);
diff = im;

%{
var = 2;
x = (-4:4);
g = exp(-x.*x/(2*var)); g  = g/sum(g);

blurred = conv2(im,g,'same');
im_b = conv2(blurred,g','same'); 

%}

for i = 1:niter
    
 % fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff;

  % North, South, East and West differences
  deltaN = diffl(1:rows,2:cols+1)   - diff;
  
  deltaS = diffl(3:rows+2,2:cols+1) - diff;
  
  deltaE = diffl(2:rows+1,3:cols+2) - diff;
  
  deltaW = diffl(2:rows+1,1:cols)   - diff;
  %deltaN = diff;deltaW;
  
  %MAD implementation
  deln_med = medfilt2(abs(deltaN),[3,3],'symmetric');
  deln_abs = abs(deltaN - deln_med);
  delN_M = 1.4826* medfilt2(deln_abs,[3,3],'symmetric') * sqrt(5);
  
%   dels_med = medfilt2(abs(deltaS),[3,3],'symmetric');
%   dels_abs = abs(deltaS - dels_med);
%   delS_M = 1.4826* medfilt2(dels_abs,[3,3],'symmetric') * sqrt(5);
    delS_M = delN_M; 
  
  delw_med = medfilt2(abs(deltaW),[3,3],'symmetric');
  delw_abs = abs(deltaW - delw_med);
  delW_M   = 1.4826* medfilt2(delw_abs,[3,3],'symmetric') * sqrt(5);
  
%   dele_med = medfilt2(abs(deltaE),[3,3],'symmetric'); 
%   dele_abs = abs(deltaE - dele_med);
%   delE_M = 1.4826* medfilt2(dele_abs,[3,3],'symmetric') * sqrt(5);
  delE_M   = delW_M;

  % Conduction

  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    
    cS = exp(-(deltaS/kappa).^2);
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
    
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
  
  elseif option == 3
      
    indN  = ( abs(deltaN) < delN_M );
    indS  = ( abs(deltaS) < delS_M );
    indW  = ( abs(deltaW) < delW_M );
    indE  = ( abs(deltaE) < delE_M );
    
    cN = zeros(size(deltaN));
    cS = zeros(size(deltaS));
    cW = zeros(size(deltaW));
    cE = zeros(size(deltaE));
    
    cN(indN) =  1/2 * (1 - ((deltaN(indN)./delN_M(indN)).^2).^2);
    cS(indS) =  1/2 * (1 - ((deltaS(indS)./delS_M(indS)).^2).^2); 
    cW(indW) =  1/2 * (1 - ((deltaW(indW)./delW_M(indW)).^2).^2);
    cE(indE) =  1/2 * (1 - ((deltaE(indE)./delE_M(indE)).^2).^2);
     
  end

  % APPLYING FOUR-POINT-TEMPLETE FOR numerical solution of DIFFUSION P.D.E.
  
  diff = diff + lambda*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW);
  
  diff_g = uint8(diff);
  
  [i, psnr(diff_g, uint8(im_good)), mse(diff_g, uint8(im_good)) ]
% figure();
% % Uncomment the following to see a progression of images
% subplot(ceil(sqrt(niter)),ceil(sqrt(niter)), i)
figure();
pause(0)
imagesc(diff), colormap(gray), axis image

end
fprintf('\n');

% Uncomment the following to see histogram of final loop coeffecient value
% figure()
% imhist(cN)
% figure()
% imhist(cS)
% figure()
% imhist(cW)
% figure()
% imhist(cE)