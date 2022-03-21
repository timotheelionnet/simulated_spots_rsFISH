% add iniconfig package to matlab path
addpath('iniconfig');

%% read baseline parameters for config file
cfg = spotGenerationParams;
cfg.buidFromFile('pairs3D.ini');

% build list of conditions 
dist = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4];
bgStd = [1,2,4];
I = 300;
sxy = [1,1.35,2];
sz = 2;

condTable = table('Size',[0,5],...
    'VariableTypes',repmat({'double'},1,5),...
    'VariableNames',{'nSpots','bgStd','I','sxy','sz'});
ctr = 1;
for i1 = 1:numel(dist)
    for i2 = 1:numel(bgStd)
        for i3 = 1:numel(I)
            for i4 = 1:numel(sxy)
                for i5 = 1:numel(sz)
                
                    fName = ['dist',num2str(dist(i1)),...
                        '_bgStd',num2str(bgStd(i2)),...
                        '_I',num2str(I(i3)),...
                        '_sxy',num2str(sxy(i4)),...
                        '_sz',num2str(sz(i5))]
                    
                    % update parameters values
                    cfg.set('distBetweenSpotPairs',[dist(i1),0]);
                    cfg.set('bgStd',bgStd(i2));
                    cfg.set('brightness',[I(i3),0]);
                    cfg.set('psf',[sxy(i4), sz(i5)]);
                    cfg.set('outFolder',['out/',fName]);
                    
                    % save parameters to config file
                    cfg.saveConfigAsIni([fName,'.ini']);
                    
                    % generate and save image series
                    cfg.generateAndSaveImageSeries(fName);
                    
                end
            end
        end
    end
end
