function img = addSpotsToImage(loc,imSize,voxSize,bgMean,bgStd,psf,...
    zMode,cutSize,mode,gain)

% function that builds a simulated image containing diffraction limited spots 


%% output
% img : a 3d or 2d image with the spots added

%% input
% imSize: array [nx ny] (2D) or [nx ny nz] (3D) setting the size of the
    % image (in pixels) in which the spots will be generated
% voxSize: physical size of a voxel vxy (2D) or [vxy,vz] (3D)
% bgMean: mean intensity value of the background
% bgStd: STDEV of the background noise (gaussian distributed)
% psf = [sigma_xy sigma_z] of the psd modeled as a 2D gaussian in physical units
% zMode: 'gaussian' means the psf is gaussian along z (recommended); 'integrated
    % gaussian' means the psf is the integral of a gaussian between planes
% cutSize:  size of the region (in units of the psf size) over which 
    % the function calculates intensity around each spot. recommended value: 3.
% mode: set to 'poisson' to model spot intensities as poisson distributed
    % (shot noise).
% gain: gain of the sensor (only used in poisson mode).

%% output
%dlocpix = [x y z I] is a nPts by 4 array with the positions of the center
%of each spot in voxel units and its integrated intensity positions  
%%
% generate gaussian distributed background stack
img = bgMean + bgStd*randn(imSize);

nDims = numel(imSize); 
% loop over all the spots; for each spot, define a bounding box surrounding
% it and compute the value of the spot psf in each pixel of the box. Then
% add the values from the box to their corresponding position in the image
% stack
for j=1 : size(loc,1)

    xmin = loc(j,1) - cutSize*psf(1)/voxSize(1); %defining the ensemble of pixels b/w xmin:xmax ymin:ymax etc
    xmax = loc(j,1) + cutSize*psf(1)/voxSize(1); %these are the pixels on which the current dot is
    ymin = loc(j,2) - cutSize*psf(1)/voxSize(1); %gonna contribute non-negligible intensity
    ymax = loc(j,2) + cutSize*psf(1)/voxSize(1); %they are defined as pixels within a 3-sigma radius from the dot
    if nDims == 3
        zmin = loc(j,3) - cutSize*psf(2)/voxSize(2);
        zmax = loc(j,3) + cutSize*psf(2)/voxSize(2);
    end
    
    xp_min = max( ceil(xmin), 1);
    xp_max = min( ceil(xmax), imSize(1) );
    yp_min = max( ceil(ymin), 1);
    yp_max = min( ceil(ymax), imSize(2) );
    if nDims == 3
        zp_min = max( ceil(zmin), 1);
        zp_max = min( ceil(zmax), imSize(3) );
    end
    
    % define a local 2D or 3D box centered on the spots
    if nDims == 3
        [xp,yp,zp] = ndgrid(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max);
    else
        [xp,yp] = ndgrid(xp_min:xp_max,yp_min:yp_max);
    end
    
    diffx1 =  (double(xp-1))*voxSize(1) - loc(j,1)*voxSize(1);
    diffx1 = diffx1 / ( sqrt(2.0)*psf(1));
    diffx2 =  (double(xp))*voxSize(1) - loc(j,1)*voxSize(1);
    diffx2 = diffx2 / ( sqrt(2.0)*psf(1));

    diffy1 =  (double(yp-1))*voxSize(1) - loc(j,2)*voxSize(1);
    diffy1 = diffy1 / ( sqrt(2.0)*psf(1) );
    diffy2 =  (double(yp))*voxSize(1) - loc(j,2)*voxSize(1);
    diffy2 = diffy2 / ( sqrt(2.0)*psf(1) );
    
    intensity = loc(j,nDims+1)* abs( erf( diffx1) - erf(diffx2) )...
                .* abs( erf( diffy1) - erf(diffy2) )/4;
            
    if nDims == 3
        if strcmp(zMode,'integrated gaussian') 
            diffz1 =  (double(zp-1))*voxSize(2) - loc(j,3)*voxSize(2);
            diffz1 = diffz1 / ( sqrt(2.0)*psf(2) );
            diffz2 =  (double(zp))*voxSize(2) - loc(j,3)*voxSize(2);
            diffz2 = diffz2 / ( sqrt(2.0)*psf(2) );

            intensity = intensity .* abs( erf( diffz1) - erf(diffz2)  );
            intensity = intensity / 2.0;
            
        elseif strcmp(zMode,'gaussian') 
            gz = exp( - ( double(zp) *voxSize(2) - loc(j,3)*voxSize(2)).^2 /(2*psf(2)^2));
            intensity = intensity .* gz;
            intensity = intensity / (psf(2)/voxSize(2)*sqrt(2*pi));
        end
    end
    
    %taking a poisson distribution if poisson mode selected
    if strcmp(mode,'poisson')
       intensity = gain*poissrnd(intensity/gain);
    end
    
    % fill in the spot intensity in the full size image
    if nDims == 3
        img(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max) = ...
            img(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max) + intensity;
    else
        img(xp_min:xp_max,yp_min:yp_max) = ...
            img(xp_min:xp_max,yp_min:yp_max) + intensity;
    end
    
end
            
end