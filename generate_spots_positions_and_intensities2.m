function loc = generate_spots_positions_and_intensities2(imSize,voxSize,cutSize,nPts,brightness,psf,generatePairs,distBetweenSpotPairs)
%% input
% imSize: array [nx ny] (2D) or [nx ny nz] (3D) setting the size of the
    % image (in pixels) in which the spots will be generated
% voxSize: physical size of a voxel vxy (2D) or [vxy,vz] (3D)
% cutSize:  defines a padding region to ensure that
          % spots aren't cropped by the boundaries of the image. If pairs
          % are generated, the padding is further extended.
% nPts: number of points to be generated
% brightness: [Itot, sdItot] the mean and sd of the integrated intensity of
    % each spot
% psf = [sigma_xy sigma_z] in physical units
% generatePairs: set to zero to generate random individual spots; set to 1
    % to generate pairs of spots with specified distance distribution
% distBetweenSpotPairs: [mean STD] array defining the mean and the std of
% the radial distances between spots in physical units (gaussian distribution of distances, random angles)

%% output
%dlocpix = [x y z I] is a nPts by 4 array with the positions of the center
%of each spot in voxel units and its integrated intensity positions  

if generatePairs == 1 % generate half the numer of points if generating pairs
    nPts = ceil(nPts/2);
    nRows = 2*nPts;
else
    nRows = nPts;
end

nDims = numel(imSize);

%% initialize the coordinates array and fill the intensity column (enforcing
% positive spot intensitiy values)
loc = zeros(nRows,nDims+1);
loc(:,nDims+1) = randn(nRows,1)*brightness(2) + brightness(1); 
loc(:,nDims+1) = max(loc(:,nDims+1),0); 
    
    
%% define the region where spot centers are allowed
% image size in physical units
lx = imSize(1)*voxSize(1); 
ly = imSize(2)*voxSize(1);
if nDims ==3
    lz = imSize(3)*voxSize(2);    
end

% padd the boundaries of the image to avoid spots getting cropped.
% in physical units
paddingSizeXY = cutSize*psf(1) + ...
    double(generatePairs)*...
    (distBetweenSpotPairs(1)+3*distBetweenSpotPairs(2));
if nDims ==3
    paddingSizeZ = cutSize*psf(2) + ...
        double(generatePairs)*...
        (distBetweenSpotPairs(1)+3*distBetweenSpotPairs(2));    
end

%% generate positions of the spot centers (or half the spot centers if creating pairs)
% in units of voxels
xRange = lx - 2*paddingSizeXY;  % in physical units
yRange = ly - 2*paddingSizeXY; 
if nDims == 3
    zRange = lz - 2*paddingSizeZ;   %max z distance of a spot from the center of the ROI or image
end

% switch to voxel units
loc(1:nPts,1) = (rand(nPts,1)*xRange + paddingSizeXY)/voxSize(1) ;
loc(1:nPts,2) = (rand(nPts,1)*yRange + paddingSizeXY)/voxSize(1) ;
if nDims == 3
    loc(1:nPts,3) = (rand(nPts,1)*zRange + paddingSizeZ)/voxSize(2) ;
end

%% add a paired spot to each spot if option is selected
if generatePairs
    % generate distribution of radial distances
    dr = distBetweenSpotPairs(2)* randn(nPts,1) + distBetweenSpotPairs(1);
    
    % generate distribution of polar and azimuth angles
    theta = pi*(rand(nPts,1)-0.5).*double(nDims==3) + pi/2; % set theta to zero if 2D
    phi = 2*pi*rand(nPts,1);
    
    % generate new positions in voxel units based on position of first spot and offset
    loc(nPts+1:nRows,1) = (voxSize(1)*loc(1:nPts,1) + dr.*sin(theta).*cos(phi))/voxSize(1);
    loc(nPts+1:nRows,2) = (voxSize(1)*loc(1:nPts,2) + dr.*sin(theta).*sin(phi))/voxSize(1);
    if nDims == 3
        loc(nPts+1:nRows,3) = (voxSize(2)*loc(1:nPts,3) + dr.*cos(theta))/voxSize(2);
    end
end
    
end