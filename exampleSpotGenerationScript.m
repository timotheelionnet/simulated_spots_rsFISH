% features some of the functionalities

% add iniconfig package to matlab path
addpath('iniconfig');

%% read parameters for config file
cfg = spotGenerationParams;
cfg.buidFromFile('configFile3D.ini'); % this is the set of parameters to simulate 3D stacks
%cfg.buidFromFile('configFile2D.ini'); % this is the set of parameters to simulate 2D images

%% build list of spot positions
% loc is an array of spot coordinates and intensities with the format
% [x y z I] (3D) or [x y I] (2D) where each row is one spot. All
% coordinates in voxel units.
% if n spot pairs are generated, spot 1 is paired with spot n+1, etc.
loc = cfg.generateListOfSpotPositions();

%% build image from list
% generate image or image stack based on parameters
img = cfg.generateSingleImageFromLoc(loc);

%% save loc file
% saved in the outFolder path defined in the parameters config file
cfg.saveLocFile(loc);

%% save parameters to config file
% saved in the outFolder path defined in the parameters config file
cfg.saveConfigAsIni('testIni.ini');

%% generate and save image series
% this is how you generate a series of simulated images
cfg.generateAndSaveImageSeries('test');

%% validate pair distances
% quick sanity check that the distances between the spot pairs follow the
% expected distribution.
if cfg.generatePairs
    nPairs = size(loc,1)/2;
    dr = loc(1:nPairs,1:3) - loc(nPairs+1:end,1:3);
    dr(:,1:2) = dr(:,1:2)*cfg.voxSize(1);
    dr(:,3) = dr(:,3)*cfg.voxSize(2);
    
    dr = sqrt(sum(dr.^2,2));
    
end