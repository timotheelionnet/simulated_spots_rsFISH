classdef spotGenerationParams < handle
    
     properties (GetAccess = 'public', SetAccess = 'public')
         outFolder; % where the images will be saved
         nImgs; % number of images to generate
         nPts; % number of spots in image
         bgMean; % mean value of the background
         bgStd; % st dev of the background noise (modeled as a gaussian distrbution)

         imSize; % size of the image to be modeled entered as [nx,ny] or [ny ny nz] for z-stack in pixel units
         voxSize; % vxy,vz (3D) vxy (2D): voxel size (in voxels)
            % set to 1 (2D) or 1,1 (3D) to only use voxel units for everything
         
         brightness; % [mean stdev] of the spot brightness distribution. 
            % spot brightnesses are gaussian distributed around the mean
         
         psf; % [sigma_xy sigma_z] of the psf modeled as a 3D gaussian 
            % (in physical units e.g. nm)

         cutSize; % size of the ROI over which the spot intensity is calculated 
            %(in units of PSF sigma)

         mode; % photon noise model. 
             % Set to '' for no noise; 
             % set to 'poisson' for poisson shot noise. When choosing poisson
             % noise, noise is modeled in each pixzel as gain*poissrnd(intensity/gain)

         gain; % gain of the sensor (only used when mode set to 'poisson').

         zMode; % how the intensity profile of the psf is modeled in z. 
             % Set to 'gaussian' for gaussian distribution (recommended); 
             % set to 'integrated gaussian' for gaussian integrated over the voxel
             % height.

         generatePairs; %set to zero to generate unique spots. 
            % Set to 1 to generate spot pairs.

         distBetweenSpotPairs; % [mean stdev] of the distribution 
            % of distances between spot pairs if generating spot pairs (in physical units e.g. nm)
         
         defaultProperties; % structure holding default values of the properties
         
         propertiesNotReadFromFile; % list of properties that shouldn't be in the config file
         
         allowedValues; % list of values allowed for some parameters
         
         iniFileName; % the path of the config file used.
         
     end
     
      methods
        %% initialize object
        function obj = spotGenerationParams()
            obj.initDefaultProperties;
            obj.initAllowedValues;
            obj.propertiesNotReadFromFile = obj.defaultProperties.propertiesNotReadFromFile;
            %obj.initProperties;
        end
        
        function initDefaultProperties(obj)
            obj.defaultProperties.outFolder = 'out';
            obj.defaultProperties.nPts = 30;
            obj.defaultProperties.nImgs = 1;
            obj.defaultProperties.bgMean = 200;
            obj.defaultProperties.bgStd = 10;
            obj.defaultProperties.imSize = [256,256,30];
            obj.defaultProperties.voxSize = [1,1];
            obj.defaultProperties.brightness = [1000,0];
            obj.defaultProperties.psf = [1.35,2];
            obj.defaultProperties.cutSize = 3;
            obj.defaultProperties.gain = 1;
            obj.defaultProperties.mode = '';
            obj.defaultProperties.zMode = 'gaussian';
            obj.defaultProperties.generatePairs = 0;
            obj.defaultProperties.distBetweenSpotPairs = 5;
            obj.defaultProperties.iniFileName = '';
            obj.defaultProperties.propertiesNotReadFromFile = ...
                {'iniFileName','propertiesNotReadFromFile',...
                'defaultProperties','allowedValues'};
            
        end
        
        function initAllowedValues(obj)
            obj.allowedValues.mode={'','poisson'};
            obj.allowedValues.zMode={'gaussian','integratedGaussian'};
            obj.allowedValues.generatePairs = [0,1];
        end
        
        function buidFromFile(obj,inputFileName)
            obj.iniFileName = inputFileName;
            if isempty(inputFileName)
                disp('config file path is empty');
                return
            end
            
            % build iniconfig object
            inConf = IniConfig();
            inConf.ReadFile(inputFileName);
            
            if isempty(inConf)
                disp(['Could not load fileSet from config file: ',...
                    obj.iniFileName]);
                status = 0;
                return
            else
                status = 1;
            end
            
            % fill in every property from the matching entry in the config
            % file
            propertyList = properties(spotGenerationParams);
            for i=1:numel(propertyList)
                if ~ismember(propertyList{i},obj.defaultProperties.propertiesNotReadFromFile)
                    p = GetValues(inConf, ...
                        'spotParameters', propertyList{i}, 'fail');
                    if ~strcmp(p,'fail') 
                        set(obj,propertyList{i},p);
                    else
                        disp(['Could not parse ',propertyList{i},' from config file, setting to default!']);
                        obj.setToDefault(propertyList{i});
                    end
                    obj.checkPropertyFormat(propertyList{i});
                end
            end
        end
        
        % function that generates an array of spot coordinates randomly generated based
        % on the parameters stored in the instance
        % the array has columns of [x y z I]. If spots are generated in
        % pairs (2*n spots total), row 1 and row n+1 belong to one pair;
        % row 2 and n+2 belong to one pair; etc.
        % loc array coordinates are in units of voxels
        function loc = generateListOfSpotPositions(obj)
            loc = generate_spots_positions_and_intensities2(...
                obj.imSize,obj.voxSize,obj.cutSize,obj.nPts,obj.brightness,...
                obj.psf,obj.generatePairs,obj.distBetweenSpotPairs);
        end
        
        % function that generates an image based on an array of spot
        % positions (loc, in [x y z I format]) and a set of parameters
        % stored in the instance properties.
        function img = generateSingleImageFromLoc(obj,loc)
            img = addSpotsToImage(loc,obj.imSize,obj.voxSize,...
                obj.bgMean,obj.bgStd,obj.psf,obj.zMode,obj.cutSize,...
                obj.mode,obj.gain);
            
        end
        
        % function that generates a series of images and saves them to outFolder 
        function generateAndSaveImageSeries(obj,fName)
            nDigits = 4;
            for i=1:obj.nImgs
                loc = obj.generateListOfSpotPositions;
                img = obj.generateSingleImageFromLoc(loc);
                curDigits = ceil(log(i+1)/log(10));
                zeroString = repmat('0',1,nDigits-curDigits);
                
                obj.saveLocFile(loc,[fName,'_',zeroString,num2str(i),'.loc']);
                obj.saveImgFile(img,[fName,'_',zeroString,num2str(i),'.tif']);
            end
        end
        
        % saves the loc file as a .loc ascii file in the designated outFolder
        % optional argument is the file name, if no name is supplied,
        % defaults to 'simulatedSpots.loc'
        function saveLocFile(obj,loc,varargin)
            if ~exist(obj.outFolder,'dir')
                mkdir(obj.outFolder);
            end
            
            if ~isempty(varargin)
                fName = varargin{1};
            else
                fName = 'simulatedSpots.loc';
            end
            save(fullfile(obj.outFolder,fName),'loc','-ascii');  
        end
        
        % saves the loc file as a .loc ascii file in the designated outFolder
        % optional argument is the file name, if no name is supplied,
        % defaults to 'simulatedSpots.loc'
        function saveImgFile(obj,img,varargin)
            if ~exist(obj.outFolder,'dir')
                mkdir(obj.outFolder);
            end
            
            if ~isempty(varargin)
                fName = varargin{1};
            else
                fName = 'simulatedSpots.tif';
            end
            save_as_tiff(img,fullfile(obj.outFolder,fName));  
        end
        
        % saves the properties of the instance as an ini config file.
        % optional argument is the name of the file (if left empty,
        % defaults to 'spotGenerationParams.ini'.
        % saved to the outFolder.
        function saveConfigAsIni(obj,varargin)
            
            if ~isempty(varargin)
                fName = varargin{1};
            else
                fName = 'spotGenerationParams.ini';
            end
            
            ini = IniConfig();
            ini.AddSections({'spotParameters'});
            propertyList = properties(spotGenerationParams);
            for i=1:numel(propertyList)
                if ~ismember(propertyList{i},...
                        obj.defaultProperties.propertiesNotReadFromFile)
                    ini.AddKeys('spotParameters', propertyList{i}, ...
                        obj.get(propertyList{i}));
                end
            end
            
            if ~exist(obj.outFolder,'dir')
                mkdir(obj.outFolder);
            end
            
            ini.WriteFile(fullfile(obj.outFolder,fName));
        end
        
        function checkPropertyFormat(obj,propertyName) 
             switch propertyName
                case 'outFolder'
                    
                case 'nPts'
                    if obj.nPts<=0
                        disp('error: nPts non positive, setting to default');
                        obj.setToDefault('nPts');
                    end
                case 'nImgs'
                    if obj.nPts<=0
                        disp('error: nImgs non positive, setting to default');
                        obj.setToDefault('nImgs');
                    end    
                case 'bgMean'
                    if obj.bgMean<=0
                        disp('error: bg Mean non positive, setting to default');
                        obj.setToDefault('bgMean');
                    end
                case 'bgStd'
                    if obj.bgStd<0
                        disp('error: bg Std non positive, setting to default');
                        obj.setToDefault('bgStd');
                    end
                case 'imSize'
                    if numel(obj.imSize)~=2 && numel(obj.imSize)~=3
                        disp('error: imSize not a 2 or 3 element array, setting to default');
                        obj.setToDefault('imSize');
                    end
                case 'voxSize'
                    if numel(obj.voxSize)~=1 && numel(obj.voxSize)~=2
                        disp('error: voxSize not a 1 or 2 element array, setting to default');
                        obj.setToDefault('voxSize');
                    end
                case 'brightness'
                    if numel(obj.brightness) == 1
                        obj.brightness  = [obj.brightness, 0];
                    end
                    
                case 'psf'
                    if numel(obj.psf)~=1 && numel(obj.psf)~=2
                        disp('error: psf not a 1 or 2 element array, setting to default');
                        obj.setToDefault('psf');
                    end
                case 'cutSize'
                    
                case 'mode'
                    if ~ismember(obj.mode,obj.allowedValues.mode)
                        disp('error: mode entry not one of the allowed flags, setting to default');
                        obj.setToDefault('mode');
                    end
                    
                case 'gain'
                    if obj.gain <=0
                        disp('error: gain is negative, setting to default');
                        obj.setToDefault('gain');
                    end
                case 'zMode'
                    if ~ismember(obj.zMode,obj.allowedValues.zMode)
                        disp('error: zMode entry not one of the allowed flags, setting to default');
                        obj.setToDefault('zMode');
                    end
                case 'generatePairs'
                    if ~ismember(obj.generatePairs,obj.allowedValues.generatePairs)
                        disp('error: generatePairs entry not one of the allowed flags, setting to default');
                        obj.setToDefault('generatePairs');
                    end
                case 'distBetweenSpotPairs'
                    
                case 'defaultProperties'
                    
                case 'propertiesNotReadFromFile'
                    
                case 'iniFileName'
                    
             end
        end
        
        function set(obj,propertyName,value) 
             switch propertyName
                case 'outFolder'
                    obj.outFolder  = value;
                case 'nPts'
                    obj.nPts  = value;
                case 'nImgs'
                    obj.nImgs  = value;
                case 'bgMean'
                    obj.bgMean  = value;
                case 'bgStd'
                    obj.bgStd  = value;
                case 'imSize'
                    obj.imSize  = value;
                case 'voxSize'
                    obj.voxSize  = value;
                case 'brightness'
                    obj.brightness  = value;
                case 'psf'
                    obj.psf  = value;
                case 'cutSize'
                    obj.cutSize  = value;
                case 'mode'
                    obj.mode  = value;
                case 'gain'
                    obj.gain  = value;
                case 'zMode'
                    obj.zMode  = value;
                case 'generatePairs'
                    obj.generatePairs  = value;
                case 'distBetweenSpotPairs'
                    obj.distBetweenSpotPairs  = value;
                case 'defaultProperties'
                    obj.defaultProperties  = value;
                case 'propertiesNotReadFromFile'
                    obj.propertiesNotReadFromFile  = value;
                case 'iniFileName'
                    obj.iniFileName  = value;
             end
        end
        
        function x = get(obj,propertyName)
            switch propertyName
                case 'outFolder'
                    x = obj.outFolder;
                case 'nPts'
                    x = obj.nPts;
                case 'nImgs'
                    x = obj.nImgs;
                case 'bgMean'
                    x = obj.bgMean;
                case 'bgStd'
                    x = obj.bgStd;
                case 'imSize'
                    x = obj.imSize;
                case 'voxSize'
                    x = obj.voxSize;
                case 'brightness'
                    x = obj.brightness;
                case 'psf'
                    x = obj.psf;
                case 'cutSize'
                    x = obj.cutSize;
                case 'mode'
                    x = obj.mode;
                case 'gain'
                    x = obj.gain;
                case 'zMode'
                    x = obj.zMode;
                case 'generatePairs'
                    x = obj.generatePairs;
                case 'distBetweenSpotPairs'
                    x = obj.distBetweenSpotPairs;
                case 'defaultProperties'
                    x = obj.defaultProperties;
                case 'propertiesNotReadFromFile'
                    x = obj.propertiesNotReadFromFile;
                case 'allowedValues'
                    x = obj.allowedValues;
                case 'iniFileName'
                    x = obj.iniFileName;
            end
        end
         
        function setToDefault(obj,propertyName)
            switch propertyName
                case 'outFolder'
                    obj.outFolder  = obj.defaultProperties.outFolder;
                case 'nPts'
                    obj.nPts  = obj.defaultProperties.nPts;
                case 'nImgs'
                    obj.nImgs  = obj.defaultProperties.nImgs;
                case 'bgMean'
                    obj.bgMean  = obj.defaultProperties.bgMean;
                case 'bgStd'
                    obj.bgStd  = obj.defaultProperties.bgStd;
                case 'imSize'
                    obj.imSize  = obj.defaultProperties.imSize;
                case 'voxSize'
                    obj.voxSize  = obj.defaultProperties.voxSize;
                case 'brightness'
                    obj.brightness  = obj.defaultProperties.brightness;
                case 'psf'
                    obj.psf  = obj.defaultProperties.psf;
                case 'cutSize'
                    obj.cutSize  = obj.defaultProperties.cutSize;
                case 'mode'
                    obj.mode  = obj.defaultProperties.mode;
                case 'gain'
                    obj.gain  = obj.defaultProperties.gain;
                case 'zMode'
                    obj.zMode  = obj.defaultProperties.zMode;
                case 'generatePairs'
                    obj.generatePairs  = obj.defaultProperties.generatePairs;
                case 'distBetweenSpotPairs'
                    obj.distBetweenSpotPairs  = obj.defaultProperties.distBetweenSpotPairs;
                case 'defaultProperties'
                    obj.defaultProperties  = obj.defaultProperties.defaultProperties;
                case 'propertiesNotReadFromFile'
                    obj.propertiesNotReadFromFile  = obj.defaultProperties.propertiesNotReadFromFile;
                case 'iniFileName'
                    obj.iniFileName  = obj.defaultProperties.iniFileName;
            end
        end
        
        function generateSetCode(obj) % lazy function to generate the text for the 'set' method code...
            propertyList = properties(spotGenerationParams);
            str = [];
            for i=1:numel(propertyList)
                curStr = ['case ''',propertyList{i},'''',newline,'    ','obj.',propertyList{i},'  = value;',newline];
                str = [str,curStr];
            end
            disp(str);
            disp(' ' );
            str = [];
            for i=1:numel(propertyList)
                curStr = ['case ''',propertyList{i},'''',newline,'    ','obj.',propertyList{i},'  = obj.defaultProperties.',propertyList{i},';',newline];
                str = [str,curStr];
            end
            disp(str);
            disp(' ');
            for i=1:numel(propertyList)
                curStr = ['case ''',propertyList{i},'''',newline,'    ','x = obj.',propertyList{i},';',newline];
                str = [str,curStr];
            end
            disp(str);
        end
        
      end
end