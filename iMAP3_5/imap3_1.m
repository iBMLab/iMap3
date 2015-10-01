function imap3_1(cfg)
%%
% iMap - Version 3.1 (March 2014)
% Free open-source Matlab toolbox for statistical mapping of eye-movement data
% http://perso.unifr.ch/roberto.caldara/index.php?page=4
% Sébastien Miellet, Yingdi Liu, Cyril R. Pernet, Guillaume A. Rousselet, Michael Papinutto, Junpeng Lao, Luca Vizioli & Roberto Caldara (2014). University of Fribourg
% Junpeng Lao (2012) for the indtorgb function. Cyril Pernet & Guillaume Rousselet (2013) for the tfce2d function
% sebastien.miellet@unifr.ch

% 1. Update details
% version 3.1
% 1- Suppression of the edge effect in the spatial smoothing of the fixation maps (Yingdi Liu).
% 2- Possibility to set the input and outpout data paths (Michael Papinutto).
% 3- No need to set the minimal number of data points required to include a pixel in the analysis. As a consequence iMap does not generate a logical mask indicating the parts of the stimulus space that will be included in the analysis.
%    Instead of using logical masks to select pixels with enough data points, we add random noise to prevent having empty cells (Yingdi Liu). 
% 4- The minimal number of unique values in the bootstrap procedures (minBootUnique) can be set individually. In the previous version it was set equal to mindatapoints (Minimal number of data points required to include a pixel in the analysis). 
% 5- iMap3.1 generates the sigShading maps showing on a single picture the significant areas corresponding to different significance threshold (Yingdi Liu).
% 6- The toolbox now generates individual fixation maps - smoothed and normalized in the stimulus space (Michael Papinutto).
% 7- All the numerical outputs are now presented in a unique txt document (iMap_report.txt).
% version 3
% 1- New "statistical engine" for iMap3. iMap3 now uses pixel-wise t-values and bootstrapped TFCE transformed scores to correct for multiple comparisons (TFCE: threshold-free cluster-enhancement; Smith & Nichols, 2009), instead of being based on the Random Field Theory.
%    The associated benefits are: a) a more appropriate and direct estimate of data variability, b) specific tests for independent and paired samples, c) a better control of the presence of false positives.
%    Many thanks to Cyril R. Pernet & Guillaume A. Rousselet for their help with the implementation of this approach.
% 2- A bug due to the impossibility of calculating effect sizes when there is no effect has been fixed.
% 3- This version includes an option to generate the raw fixation maps.
% 4- This version does no longer create temporary working files.
% 5- A normalization of the individual maps is done across the search space in order to represent the individual fixation bias.
% 6- Improved appearance of the statistical fixation maps.
% 7- Possibility to set the transparency of the background image and the statistical maps.
%
% 2. Input data
% An input data file is a matrix with a fixation per line. The only crucial data are the coordinates (x and y) and duration of the fixations and the item numbers. Any other column can be used for specifying your conditions.
% Typically the data files are .mat files (called data1.mat, data2.mat,...) containing matrices called "summary". This can be obtained from any txt file (e.g. fixation report from EyeLink, Data Viewer, custom preprocessing code,...)
% See the examples for the input data format.
%
% 3. How to use iMap3_1
% Copy imap3_1 and its support functions in the folder containing the data file (or alternatively set paths for input and output data).
% Set the parameters (see configuration structure section) in a configuration structure (e.g. cfg.xSize=400). Default values will be used for non-specified parameters.
% Call the imap3_1 function with the configuration structure (e.g. imap3_1(cfg)).
% See examples at the end of this document and the example folders for more specific explanations.
% Please keep in mind that running iMap takes a bit of time due to the use of bootstrapping and TFCE. We included wait bars so the user can keep track of the analysis progression.
% Some analyses with small stimulus size (or downsampled stimulus space)can be performed on a 8Go RAM computer. However, we strongly recommend using a 16 or even better 32Go RAM computer. This code was tested on a 72Go node.
% We recommend using iMap3_1 at first with the default settings to have a general view on the results. Then, in a second step for the final analysis, we advise to use a higher number of bootstraps (1000 or more for better estimate) and, if necessary, to set the color scale and transparency of the maps.
%
% 4. Configuration structure
% VARIABLES that can be set in the cfg structure
% e.g. cfg.xSize=400. See examples at the end of this document.
%
% 1-	xSize and ySize: stimulus size in pixels (e.g. 382, 390)
% IMPORTANT: Please keep in mind that the stimuli dimensions (xSize and ySize) might be inverted depending on whether the user considers them in graph coordinates (abscissa/ordinate, bottom left origin), screen coordinates (top left origin) or matrices (number of lines first, number of columns second). Here we consider matrix coordinates.
% 2-	columnx, columny, columnduration, columnitem: specify the column number for x, y coordinates, fixation durations and item number. This allows flexible data format. By defaults these columns are 1, 2, 3 and 4
% 3-	dataset1 and dataset2: specify the data .mat files that will be tested/compared. For example [1:20], [21:40] to compare data1 to data20 with data 21 to data40. The second data set is optional. If only one dataset is tested, iMap produces a statistical map and eye-tracking indexes for this dataset. If two datasets are specified, iMap provides the statistical maps and eye-tracking indexes for both dataset and the difference map and indexes.
% 4-	twosampletest: 1=paired or 2=independent
% 5-	smoothingpic: Standard deviation in pixels of the Gaussian used for the data map smoothing. The default value is 10 pixels.
% IMPORTANT: Please note that the smoothing should correspond to the actual viewing conditions (resolution, size, distance of the screen) and spatial resolution of the eye tracker.
% 6-	maptype: 1 for fixation duration maps, 2 for number of fixations maps. The default value is 1.
% 7-	firstfix: This option allows to ignore the first fixation of each trial. This is particularly useful if the stimuli are centred and a central fixation cross is presented before the trials. 1 (default option) keeps all the fixations, 2 ignores the first fixation of each trial.
% 8-	backgroundfile: e.g. 'facebackground.tif'. This option allows adding a background picture to the statistical fixation maps.
% 9-	specificfix: To select one or several specific fixations. e.g. [3 3] or [1 3]. This value is optional.
% 10-	searchspace: By default the stimulus space, xSize * ySize.  The search space size can be specified with a logical mask, i.e. by entering the file name of a black and white picture (e.g. "facemask.tif") where the white part indicates the search space.
% 11-	scaledownup: To be specified as a 2 values vector ([scaledown scaleup]). It allows to set the color coded scale. It does present the advantage to have the same scale for the individual (specific to datasets) maps and the contrast map. We recommend to run iMap a first time without setting this parameter in order to get an idea of the range of the t-values and then to set it in a second time.
% 12-	sigthres: Significativity threshold (for pixel-test and cluster-test). One-tailed for the individual maps, two-tailed for the difference map. By default .05 (.025 for the contrast map)
% 13-	nboot: Number of resamples for the multiple comparisons correction (default 500).
% 14-	rawmaps: Generates tiff images (called rawfix1 and rawfix2) of the raw (without normalization and smoothing) fixation locations and durations. 1 = no, 2 = yes (default)
% 15-	transpim: Set the transparency of the background image from 0 to 1 (default = 1)
% 16-	transpmap: Set the transparency of the statistical map from 0 to 1 (default = .65)
% 17-   minBootUnique: Minimal number of unique values in the bootstrap procedures (default = 2)
% 18-   avgnum: Number of iterations for random noise maps that will be averaged (default = 50)
% 19-   pvalues: Vector containing a series of significance thresholds that will be used to generate the sigShading maps (default = [0.005 0.01 0.025 0.05 0.1 0.2]) 
%
% 5. Output
% 1-  Fixation maps: iMap3 creates .tif pictures of the single and difference fixation maps merged with a background picture. It displays the significant areas (displayed heat maps) after multiple comparisons correction. The color coding of the heat maps indicates the t-values. It also has the options to create .tif pictures with normalized scales and raw (without smoothing and normalization) fixation maps. Please see examples.
% 2-  Raw fixation maps: Raw (without normalization and smoothing) fixation locations and durations for each data set, averaged across participants. The raw maps are generated by default.
% 3-  Individual fixation maps: individual fixation maps (smoothed and normalized in the stimulus space) located in a folder called Individualmaps
% 4-  sigShading maps showing on a single picture the significant areas corresponding to different significance threshold.
% 5-  All the numerical outputs are now presented in a unique txt document
% (iMap_report.txt). These numerical outputs include:
%     a-  Global eye-tracking measures for each participant: number of fixations, total fixation duration, mean fixation duration, path length and mean saccade length. 
%     b-  Eye-tracking measures in significant areas (and the rest of the picture): mean fixation duration for area 1 then for area 2 then for the rest of the picture. In the same way are path length, total fixation duration and number of fixations. Please refer to the code for the exact output structure that might vary depending on the observed significant effects.
%     c-  Z-scores: The iMap toolbox creates a text file called Zscore.txt that include the mean Zscores in the significant area for (respective columns) the dataset 1 in the area 1 and area 2 (areas in which the fixation durations are significantly longer for dataset 1 and 2 respectively), dataset 2 in the area 1 and area 2. Please refer to the code for the exact output structure that might vary depending on the observed significant effects.
%     d-  Effect sizes: It also creates a .txt file with the cohen's d between both datasets for area 1 and 2 (areas in which the fixation durations are significantly longer for dataset 1 and 2 respectively). The file is called cohend.txt. The same data set is used to both determine the significant areas with iMap and compute the effect sizes from these significant areas (for selection and selective analysis). Henceforth, the non-independent selective analysis might distort descriptive statistics (see Kriegeskorte, Simmons, Bellgowan & Baker, 2009; example 2). It is relatively straightforward to generate the distribution of independent random split-data effect sizes computed with a resampling procedure. This option is not implemented by default in iMap as it requires a very long processing time.
%   
% 8. Credits
% tfce2d is adapted from LIMO-EEG.
%    Pernet, C.R., Chauveau, N., Gaspar, C.M., & Rousselet, G.G. (2011). LIMO EEG: A Toolbox for Hierarchical LInear MOdeling of ElectroEncephaloGraphic Data. Computational Intelligence and Neuroscience. Article ID 831409, doi:10.1155/2011/831409
%
% 9. Disclaimer
% iMap is a free software; you can redistribute it and/or modify it.
% We cannot be hold responsible for any damage that may (appear to) be caused by the use of iMap. Use at your own risk.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% Please cite us and let us know if you have any comment or suggestion.
% Thank you to the users who sent us feedbacks.

% CITATION
% If you use iMap, please cite us. You could use the following sentence:
% Eye movement data were analysed with iMap3.1 (Caldara and Miellet, 2011), which implements TFCE for signal enhancement (Pernet, Chauveau, Gaspar, Rousselet, 2011; Smith and Nichols, 2009).
%
% References:
% Caldara, R., & Miellet, S. (2011). iMap: A Novel Method for Statistical Fixation Mapping of Eye Movement data. Behavior Research Methods, 43(3), 864-78
% Pernet, C.R., Chauveau, N., Gaspar, C.M., & Rousselet, G.G. (2011). LIMO EEG: A Toolbox for Hierarchical LInear MOdeling of ElectroEncephaloGraphic Data. Computational Intelligence and Neuroscience. Article ID 831409, doi:10.1155/2011/831409
% Smith S.M., & Nichols, T.E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), 83-98. doi: 10.1016/j.neuroimage.2008.03.061.

%% Reading the cfg structure and setting default values

noiseLevel=0.5; 

if isfield(cfg,'inputPath')
    inputPath=cfg.inputPath;
else
    inputPath=pwd;
end
if isfield(cfg,'outputPath')
    outputPath=cfg.outputPath;
else
    outputPath=pwd;
end
if isfield(cfg, 'pvalues')
    pvalues=cfg.pvalues;
else
    pvalues=[0.005 0.01 0.025 0.05 0.1 0.2];
end
if isfield(cfg, 'minBootUnique')
    minBootUnique=cfg.minBootUnique;
else
    minBootUnique=2;
end
if isfield(cfg, 'avgnum')
    avgnum=cfg.avgnum;
else
    avgnum=50;
end
if isfield(cfg,'xSize')
    xSize=cfg.xSize;
else
    error('imap needs the size of the stimulus. Please specify cfg.xSize');
end
if isfield(cfg,'ySize')
    ySize=cfg.ySize;
else
    error('imap needs the size of the stimulus. Please specify cfg.ySize');
end
if isfield(cfg,'columnx')
    columnx=cfg.columnx;
else
    columnx=1;
end
if isfield(cfg,'columny')
    columny=cfg.columny;
else
    columny=2;
end
if isfield(cfg,'columnduration')
    columnduration=cfg.columnduration;
else
    columnduration=3;
end
if isfield(cfg,'columnitem')
    columnitem=cfg.columnitem;
else
    columnitem=4;
end
if isfield(cfg,'dataset1')
    dataset1=cfg.dataset1;
else
    error('imap needs at least 1 dataset. Please specify cfg.dataset1');
end
if isfield(cfg,'dataset2')
    numberofdataset=2;
    dataset2=cfg.dataset2;
    if isfield(cfg,'twosampletest') % twosampletest: 1=paired or 2=independant
        twosampletest=cfg.twosampletest;
    else
        error('please indicate if the samples are paired(1) or independant(2) in cfg.twosampletest');
    end
else
    numberofdataset=1;
    dataset2=[];
end
if isfield(cfg,'specificfix')
    firstfix=cfg.specificfix(1);
else
    firstfix=1;
end
if isfield(cfg,'smoothingpic')
    smoothingpic=cfg.smoothingpic;
else
    smoothingpic=10;
end
if isfield(cfg,'maptype')
    maptype=cfg.maptype;
else
    maptype=1;
end
if isfield(cfg,'firstfix')
    firstfix=cfg.firstfix;
else
    firstfix=1;
end
if isfield(cfg,'backgroundfile')
    backgroundfile=cfg.backgroundfile;
else
    backgroundfile=[];
end
if isfield(cfg,'sigthres')
    p=cfg.sigthres;
else
    p=.05;
end
if isfield(cfg,'nboot')
    nboot=cfg.nboot;
else
    nboot=500;
end
if isfield(cfg,'rawmaps') % 1 = no, 2 = yes (default)
    rawmaps=cfg.rawmaps;
else
    rawmaps=2;
end
if isfield(cfg,'transpim')
    transpim=cfg.transpim;
else
    transpim=1;
end
if isfield(cfg,'transpmap')
    transpmap=cfg.transpmap;
else
    transpmap=.65;
end

datatotal=[dataset1 dataset2];

imapreport = fopen('iMap_report.txt','a');
fprintf(imapreport,'%11s %12s\n','iMap report',date);
fprintf(imapreport,'%64s\n \n \n','----------------------------------------------------------------');

%% Global eye-tracking measures for each dataset
fprintf(imapreport,'%-64s\n','Global eye-tracking measure (in full stimulus space)');
for datasetnb=1:numberofdataset
    fprintf(imapreport,'%-8s\n','      ');
    if datasetnb==1
        dataset=dataset1;
        fprintf(imapreport,'%-8s\n','Dataset1');
    elseif datasetnb==2
        dataset=dataset2;
        fprintf(imapreport,'%-8s\n','Dataset2');
    end
    nbfixtrial=[];
    totalfixdur=[];
    meanfixdur=[];
    pathlength=[];
    meansacclength=[];
    for datafilenumber=1:length(dataset)
        summary=[];
        datatoload=[inputPath '\data' num2str(dataset(datafilenumber))];
        load(datatoload); % The name of the matrix is 'summary'
        [nbfix, ~]=size(summary);
        cumulfixdur=[];
        cumulsaccadelength=[];
        numfix=[];
        cumulfixdur(1)=summary(1,columnduration);
        cumulsaccadelength(1)=0;
        numfix(1)=1;
        trialnb=0;
        for fix=2:nbfix
            if summary(fix,columnitem)==summary(fix-1,columnitem)
                numfix(fix)=numfix(fix-1)+1;
                cumulfixdur(fix)=cumulfixdur(fix-1)+summary(fix,columnduration);
                cumulsaccadelength(fix)=cumulsaccadelength(fix-1)+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
            elseif summary(fix,columnitem)~=summary(fix-1,columnitem)
                trialnb=trialnb+1;
                nbfixtrial(trialnb,datafilenumber)=numfix(fix-1);
                totalfixdur(trialnb,datafilenumber)=cumulfixdur(fix-1);
                meanfixdur(trialnb,datafilenumber)=cumulfixdur(fix-1)/numfix(fix-1);
                pathlength(trialnb,datafilenumber)=cumulsaccadelength(fix-1);
                meansacclength(trialnb,datafilenumber)=cumulsaccadelength(fix-1)/numfix(fix-1);
                cumulfixdur(fix)=summary(fix,columnduration);
                cumulsaccadelength(fix)=0;
                numfix(fix)=1;
            end
        end
        trialnb=trialnb+1;
        nbfixtrial(trialnb,datafilenumber)=numfix(fix);
        nbfixtrial(nbfixtrial==0)=NaN;
        totalfixdur(trialnb,datafilenumber)=cumulfixdur(fix);
        totalfixdur(totalfixdur==0)=NaN;
        meanfixdur(trialnb,datafilenumber)=cumulfixdur(fix)/numfix(fix);
        meanfixdur(meanfixdur==0)=NaN;
        pathlength(trialnb,datafilenumber)=cumulsaccadelength(fix);
        pathlength(pathlength==0)=NaN;
        meansacclength(trialnb,datafilenumber)=cumulsaccadelength(fix)/numfix(fix);
        meansacclength(meansacclength==0)=NaN;
    end
    nbfixsbj=(nanmean(nbfixtrial))';
    totalfixdursbj=(nanmean(totalfixdur))';
    meanfixdursbj=(nanmean(meanfixdur))';
    pathlengthsbj=(nanmean(pathlength))';
    meansacclengthsbj=(nanmean(meansacclength))';
       
    eyebasicdataset=[nbfixsbj totalfixdursbj.*1000 meanfixdursbj.*1000 pathlengthsbj meansacclengthsbj];
    fprintf(imapreport,'%5s %12s %11s %14s %13s\n','NbFix','TotalFix(ms)','MeanFix(ms)','TotalPath(pix)','MeanSacc(pix)');
    fprintf(imapreport,'%5.2f %12.0f %11.0f %14.0f %13.0f\n',eyebasicdataset');
    nbfixsbj=[];
    totalfixdursbj=[];
    meanfixdursbj=[];
    pathlengthsbj=[];
    meansacclengthsbj=[];
    if datasetnb==numberofdataset
        fprintf(imapreport,'%64s\n \n \n','----------------------------------------------------------------');
    end
end
clear ans summary datatoload nbfixsbj pathlengthsbj cumulfixdur meansacclength cumulsaccadelength dataset eyebasicdataset datasetnb datafilenumber meanfixdur meanfixdursbj meansacclentgh meansacclengthsbj nbfix nbfixtrial nbvariables numfix pathlength totalfixdur totalfixdursbj trialnb

%% Descriptive fixation maps
% Original fixation duration matrix (before adding noise)
origmatrixdurationtotal=zeros(length(datatotal),xSize, ySize);
noiseParams=zeros(length(datatotal),1);

for datafilenumber=1:length(datatotal)
    summary=[];
    % POINTS step: Creating a matrix x by y (stimulus size) with the cumulated fixation durations for each pixel
    datatoload=[inputPath '\data' num2str(datatotal(datafilenumber))];
    load(datatoload); % The name of the matrix is 'summary'
    [nbfix nbvariables]=size(summary);
    if isfield(cfg,'specificfix')
        nbfix=cfg.specificfix(2);
    end
    matrixduration = zeros(xSize, ySize);
    for fix=firstfix:nbfix
        if firstfix>=2
            if summary(fix, columnitem)==summary(fix-1, columnitem)
                coordX = round(summary(fix, columny)); % Here we swap x and y (difference between screen coordinates and matrix [first number=lines, second=columns])
                coordY = round(summary(fix, columnx));
                if coordX<xSize && coordY<ySize && coordX>0 && coordY>0 % Here we consider only the fixations inside the stimulus space
                    if maptype==1
                        matrixduration(coordX, coordY) = matrixduration(coordX, coordY) + summary(fix, columnduration);
                    elseif maptype==2
                        matrixduration(coordX, coordY) = matrixduration(coordX, coordY) + 1;
                    end
                end
            end
        elseif firstfix==1
            coordX = round(summary(fix, columny));
            coordY = round(summary(fix, columnx));
            if coordX<xSize && coordY<ySize && coordX>0 && coordY>0
                if maptype==1
                    matrixduration(coordX, coordY) = matrixduration(coordX, coordY) + summary(fix, columnduration);
                elseif maptype==2
                    matrixduration(coordX, coordY) = matrixduration(coordX, coordY) + 1;
                end
            end
        end
    end
    origmatrixdurationtotal(datafilenumber,:,:)=matrixduration;
    % Individual noise level
    nonzero=sort(matrixduration(matrixduration~=0))';
    noiseParams(datafilenumber)=mean(nonzero(1:round(length(nonzero)*0.05)))*noiseLevel;
end

% Generate raw fixation map pictures (averaged across observers)
if rawmaps==2
    for datasetnb=1:numberofdataset
        if datasetnb==1
            dataset=dataset1;
            deb=0;
        elseif datasetnb==2
            dataset=dataset2;
            deb=length(dataset1);
        end
        poolraw = zeros(xSize, ySize);
        for ii=1:length(dataset)
            poolraw=squeeze(origmatrixdurationtotal(ii+deb,:,:))+poolraw;
        end
        figure; imagesc(poolraw);
        if datasetnb==1
            nametosave = [outputPath '\rawfix1'];
        elseif datasetnb==2
            nametosave = [outputPath '\rawfix2'];
        end
        print('-dtiff ','-r72',nametosave)
        close(gcf);
    end
end
clear  filtered_mat nametosave f_mat f_fil x y smoothpic gaussienne cumulfixdur cumulsaccadelength matrixduration summary poolraw deb coordX coordY datatoload datafilenumber datasetnb fix ii nbfix tmp;

%% Individual maps (smoothed and normalized, without noise), stored in specific subfolders
matrixdurationtotal=origmatrixdurationtotal;
smoothpictotal=zeros(size(matrixdurationtotal));
if mod(xSize,2)==0 && mod(ySize,2)==0
    [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5);
elseif mod(xSize,2)==1 && mod(ySize,2)==0
    [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5+1);
elseif mod(xSize,2)==0 && mod(ySize,2)==1
    [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5);
elseif mod(xSize,2)==1 && mod(ySize,2)==1
    [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5+1);
end
gaussienne = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
gaussienne = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
f_fil = fft2(gaussienne);

for ii=1:size(matrixdurationtotal,1)
    A=squeeze(matrixdurationtotal(ii,:,:));
    B=flipdim(A,1);
    C=flipdim(A,2);
    D=flipdim(C,1);
    total=[A C; B D];
    f_mat = fft2(total); % 2D fourrier transform on the Symmetrized points matrix
    filtered_mat = f_mat .* f_fil;
    smoothpictmp=real(fftshift(ifft2(filtered_mat)));
    smoothpic = smoothpictmp(1:xSize, 1:ySize);
    tmp=(smoothpic - mean(smoothpic(:)))/std(smoothpic(:));
    tmp=tmp-min(tmp(:));
    tmp(tmp<0.05)=0;
    smoothpictotal(ii,:,:)=tmp;
end
nametosave = [outputPath '\Individualmaps'];
mkdir(nametosave);
nametosave = [outputPath '\Individualmaps\Indmapsdataset1'];
mkdir(nametosave);

for ii = 1:length(dataset1)
    Indivmap=squeeze(smoothpictotal(ii,:,:));
    h=figure('Visible','off');
    imagesc(Indivmap); colorbar;
    nametosave = [outputPath '\Individualmaps\Indmapsdataset1\Indmappart' num2str(ii) 'set1.tiff'];
    saveas(h,nametosave,'tif');
    clear Indivmap h; close all;
end

if numberofdataset==2
    nametosave = [outputPath '\Individualmaps\Indmapsdataset2'];
    mkdir(nametosave);
    jj=0;
    for ii=length(dataset1)+1:length(datatotal)
        jj=jj+1;
        Indivmap=squeeze(smoothpictotal(ii,:,:));
        h=figure('Visible','off');
        imagesc(Indivmap); colorbar;
        nametosave = [outputPath '\Individualmaps\Indmapsdataset2\Indmappart' num2str(jj) 'set2.tiff'];
        saveas(h,nametosave,'tif');
        clear Indivmap h; close all;
    end
end
clear A B C D nametosave ii jj datafilenumber f_fil f_mat filtered_mat gaussienne matrixduration smoothpic smoothpictmp tmp total x xSize2 y ySize2  smoothpictotal Indivmap
close all;

%% Generating observed t-maps: Adding uniformly distributed random noise (avgnum times), averaging the resulting t-maps and applying TFCE
% We add uniformly distributed random noise on each individual datafile and compute t-maps
% We repeat this procedure n=avgnum times before averaging (thus cancelling out the random noise) and applying TFCE
obstmapdataset1total=zeros(avgnum,xSize,ySize);
obstmapdataset2total=zeros(avgnum,xSize,ySize);
obstmapcontrasttotal=zeros(avgnum,xSize,ySize);

rng('shuffle'); % seeding on clock
hwaitbar=waitbar(0,'Averaging TFCEd t-maps');
for ii=1:avgnum
    waitbar(ii / avgnum)   
    % add noise indepdently for each datafile
    noise=zeros(length(datatotal),xSize,ySize);
    for datafile=1:size(origmatrixdurationtotal,1)
        noise(datafile,:,:)=rand(xSize,ySize)*noiseParams(datafile);
    end
    matrixdurationtotal=origmatrixdurationtotal+noise; 
    % smoothing
    smoothpictotal=zeros(size(matrixdurationtotal));
    if mod(xSize,2)==0 && mod(ySize,2)==0
        [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5);
    elseif mod(xSize,2)==1 && mod(ySize,2)==0
        [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5+1);
    elseif mod(xSize,2)==0 && mod(ySize,2)==1
        [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5);
    elseif mod(xSize,2)==1 && mod(ySize,2)==1
        [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5+1);
    end
    gaussienne = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
    gaussienne = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
    f_fil = fft2(gaussienne);
    for jj=1:size(matrixdurationtotal,1)
        A=squeeze(matrixdurationtotal(jj,:,:));
        B=flipdim(A,1);
        C=flipdim(A,2);
        D=flipdim(C,1);
        total=[A C; B D];
        f_mat = fft2(total); % 2D fourrier transform on the Symmetrized points matrix
        filtered_mat = f_mat .* f_fil;
        smoothpictmp=real(fftshift(ifft2(filtered_mat)));
        smoothpic = smoothpictmp(1:xSize, 1:ySize);
        tmp=(smoothpic - mean(smoothpic(:)))/std(smoothpic(:));
        smoothpictotal(jj,:,:)=tmp;
    end  
    % t-map for dataset1
    [~,~,~,stats]=ttest((smoothpictotal(1:length(dataset1),:,:)),mean(mean(mean((smoothpictotal(1:length(dataset1),:,:))))));
    dataset1tmp=squeeze(stats.tstat);
    dataset1tmp(dataset1tmp<0)=0;
    obstmapdataset1total(ii,:,:)=dataset1tmp;   
    % t-map for dataset2
    if ~isempty(dataset2)
        [~,~,~,stats]=ttest((smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:)),mean(mean(mean((smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:))))));
        dataset2tmp=squeeze(stats.tstat);
        dataset2tmp(dataset2tmp<0)=0;
        obstmapdataset2total(ii,:,:)=dataset2tmp;    
        % t-map for contrast 
        if twosampletest==1
            [~,~,~,stats]=ttest(smoothpictotal(1:length(dataset1),:,:),smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:));
            obstmapcontrasttotal(ii,:,:)=squeeze(stats.tstat);
        else
            [~,~,~,stats]=ttest2((smoothpictotal(1:length(dataset1),:,:)),(smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:)),p,'both','unequal',1);
            obstmapcontrasttotal(ii,:,:)=squeeze(stats.tstat);
        end
    end
end

% averaged and TFCEd t-maps
obstmapdataset1=squeeze(mean(obstmapdataset1total,1));
obstfcedataset1=tfce2d(obstmapdataset1);
if ~isempty(dataset2)
    obstmapdataset2=squeeze(mean(obstmapdataset2total,1));
    obstfcedataset2=tfce2d(obstmapdataset2);
    obstmapcontrast=squeeze(mean(obstmapcontrasttotal,1));
    obstmapcontrastpos=obstmapcontrast.*(obstmapcontrast>0);
    obstmapcontrastneg=obstmapcontrast.*(obstmapcontrast<0);
    obstfcecontrastpos=tfce2d(obstmapcontrastpos);
    obstfcecontrastneg=tfce2d(-obstmapcontrastneg);
end
close(hwaitbar)
clear obstmapdataset1total obstmapdataset2total obstmapcontrasttotal
% clear smoothpictotal gaussienne f_fil f_mat filtered_mat

%% Bootstrapping for multiple comparisons correction
bootAvgDatset1Total=zeros(nboot,xSize,ySize);
bootAvgDatset2Total=zeros(nboot,xSize,ySize);
bootAvgContrastPosTotal=zeros(nboot,xSize,ySize);
bootAvgContrastNegTotal=zeros(nboot,xSize,ySize);

boot_table1 = randi(length(dataset1),length(dataset1),nboot*3);
if ~isempty(dataset2)
    if twosampletest==2
        boot_table2 = randi(length(dataset2),length(dataset2),nboot*3);
    else
        boot_table2 = boot_table1;
    end
end
B1=1; B2=1;

hwaitbar=waitbar(0,'Bootstrapping...');
for bb=1:nboot
    waitbar(bb/nboot)
    % for dataset1
    while length(unique(boot_table1(:,B1))) < minBootUnique
        B1=B1+1;
    end
    % for dataset2
    while length(unique(boot_table2(:,B2))) < minBootUnique
        B2=B2+1;
    end
    
    actualBoots(1,bb)=B1; 
    actualBoots(2,bb)=B2; 
    
    % prelocate matrices for to-be-averaged t-maps
    avgDataset1tmaps=zeros(avgnum,xSize,ySize);
    avgDataset2tmaps=zeros(avgnum,xSize,ySize);
    avgContrastmapsPos=zeros(avgnum,xSize,ySize);
    avgContrastmapsNeg=zeros(avgnum,xSize,ySize);
    
    
    % add noise for 'avgnum' times
    for avg=1:avgnum
        
        noise=zeros(length(datatotal),xSize,ySize);
        for datafile=1:size(origmatrixdurationtotal,1)
            noise(datafile,:,:)=rand(xSize,ySize)*noiseParams(datafile);
        end
        matrixdurationtotal=origmatrixdurationtotal+noise;
%         matrixdurationtotal=origmatrixdurationtotal+rand(size(origmatrixdurationtotal))*noiseParam;
        

        smoothpictotal=zeros(size(matrixdurationtotal));
        if mod(xSize,2)==0 && mod(ySize,2)==0
            [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5);
        elseif mod(xSize,2)==1 && mod(ySize,2)==0
            [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5, -floor(xSize)+.5:floor(xSize)-.5+1);
        elseif mod(xSize,2)==0 && mod(ySize,2)==1
            [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5);
        elseif mod(xSize,2)==1 && mod(ySize,2)==1
            [x, y] = meshgrid(-floor(ySize)+.5:floor(ySize)-.5+1, -floor(xSize)+.5:floor(xSize)-.5+1);
        end
        gaussienne = exp(- (x .^2 / smoothingpic ^2) - (y .^2 / smoothingpic ^2));
        gaussienne = (gaussienne - min(gaussienne(:))) / (max(gaussienne(:)) - min(gaussienne(:)));
        f_fil = fft2(gaussienne);
        for ii=1:size(matrixdurationtotal,1)
            A=squeeze(matrixdurationtotal(ii,:,:));
            B=flipdim(A,1);
            C=flipdim(A,2);
            D=flipdim(C,1);
            total=[A C; B D];
            f_mat = fft2(total); % 2D fourrier transform on the Symmetrized points matrix
            filtered_mat = f_mat .* f_fil;
            smoothpictmp=real(fftshift(ifft2(filtered_mat)));
            smoothpic = smoothpictmp(1:xSize, 1:ySize);
            tmp=(smoothpic - mean(smoothpic(:)))/std(smoothpic(:));
            tmp=tmp-min(tmp(:));
            tmp(tmp<0.05)=0;
            smoothpictotal(ii,:,:)=tmp;
        end
        
        clear A B C D total f_fil f_mat filtered_mat smoothpictmp tmp x y
        
        % t-test on the smoothed-noised-images
        centeredgp1=(smoothpictotal(1:length(dataset1),:,:))-repmat(mean((smoothpictotal(1:length(dataset1),:,:)),1),[length(dataset1) 1 1]);
        dataset1map=squeeze(mean((smoothpictotal(1:length(dataset1),:,:)),1));
        
        [~,~,~,stats]=ttest(centeredgp1(boot_table1(:,B1),:,:),p);
        dataset1tmp=squeeze(stats.tstat);
        dataset1tmp(dataset1tmp<0)=0;
        avgDataset1tmaps(avg,:,:)=dataset1tmp;
        
        if ~isempty(dataset2)
            centeredgp2=(smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:))-repmat(mean((smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:)),1),[length(dataset2) 1 1]);
            dataset2map=squeeze(mean((smoothpictotal(length(dataset1)+1:length(dataset1)+length(dataset2),:,:)),1));
            
            [~,~,~,stats]=ttest(centeredgp2(boot_table2(:,B2),:,:),p);
            dataset2tmp=squeeze(stats.tstat);
            dataset2tmp(dataset2tmp<0)=0;
            avgDataset2tmaps(avg,:,:)=dataset2tmp;
            
            if twosampletest==1 % twosampletest: 1=paired
                [~,~,~,stats]=ttest(centeredgp1(boot_table1(:,B1),:,:),centeredgp2(boot_table2(:,B2),:,:));
            else  % 2=independent
                [~,~,~,stats]=ttest2(centeredgp1(boot_table1(:,B1),:,:),centeredgp2(boot_table2(:,B2),:,:),p,'both','unequal',1);
            end
            avgContrastTmap=squeeze(stats.tstat);
            avgContrastmapsPos(avg,:,:)=avgContrastTmap.*(avgContrastTmap>0);
            avgContrastmapsNeg(avg,:,:)=avgContrastTmap.*(avgContrastTmap<0);
        end
    end
    
    clear matrixdurationtotal
    
    bootAvgDatset1Total(bb,:,:)=squeeze(mean(avgDataset1tmaps,1));
    bootAvgDatset2Total(bb,:,:)=squeeze(mean(avgDataset2tmaps,1));
    if ~isempty(dataset2)
        bootAvgContrastPosTotal(bb,:,:)=squeeze(mean(avgContrastmapsPos,1));
        bootAvgContrastNegTotal(bb,:,:)=squeeze(mean(avgContrastmapsNeg,1));
    end
    
    clear avgDataset1tmaps avgDataset2tmaps avgContrastmapsPos avgContrastmapsNeg
    
    B1=B1+1;
    B2=B2+1;
    
    display(['boot=', num2str(bb), ' finished'])
end
close(hwaitbar)
clear origmatrixdurationtotal

bootTfceDataset1=tfce2d(bootAvgDatset1Total);
% clear bootAvgDatset1Total
bootTfceDataset2=tfce2d(bootAvgDatset2Total);
% clear bootAvgDatset2Total

if ~isempty(dataset2)
    bootTfceContrastPos=tfce2d(bootAvgContrastPosTotal);
    bootTfceContrastNeg=tfce2d(bootAvgContrastNegTotal);
end
% clear  bootAvgContrastPosTotal bootAvgContrastNegTotal


%% sorting TFCE scores
max_tfceH0_dataset1=zeros(nboot,1);
max_tfceH0_dataset2=zeros(nboot,1);
max_tfceH0_pos=zeros(nboot,1);
max_tfceH0_neg=zeros(nboot,1);
for b=1:nboot
    tmp1=squeeze(bootTfceDataset1(b,:,:));
    tmp2=squeeze(bootTfceDataset2(b,:,:));
    max_tfceH0_dataset1(b) = max(tmp1(:));
    max_tfceH0_dataset2(b) = max(tmp2(:));
    if ~isempty(dataset2)
        tmpPos=squeeze(bootTfceContrastPos(b,:,:));
        tmpNeg=squeeze(bootTfceContrastNeg(b,:,:));
        max_tfceH0_pos(b) = max(tmpPos(:));
        max_tfceH0_neg(b) = max(tmpNeg(:));
    end
end
clear tmp1 tmp2 tmpPos tmpNeg

v1=sort(max_tfceH0_dataset1);
v2=sort(max_tfceH0_dataset2);
if ~isempty(dataset2)
    vpos=sort(max_tfceH0_pos);
    vneg=sort(max_tfceH0_neg);
end

figure;
subplot(2,1,1); plot(v1);
subplot(2,1,2); plot(v2);
figure;
subplot(2,1,1); plot(vpos);
subplot(2,1,2); plot(vneg);



%% shading significance map(s)

% shading significance map of dataset1
all_tfce_sig_dataset1=zeros(xSize,ySize,length(pvalues));
for ii=1:length(pvalues)
    thisp=pvalues(ii);
    all_tfce_sig_dataset1(:,:,ii)= obstfcedataset1> v1(round((1-thisp)*nboot));
end

sigShades_dataset1=zeros(xSize,ySize,length(pvalues));
sigShades_dataset1(:,:,1)=all_tfce_sig_dataset1(:,:,1);
for ii=2:length(pvalues)
    sigShades_dataset1(:,:,ii)=(ones(xSize,ySize)-all_tfce_sig_dataset1(:,:,ii-1)).*all_tfce_sig_dataset1(:,:,ii)*ii+sigShades_dataset1(:,:,ii-1);
end

figure;
% colormap
map=jet;
newmap=map(1:floor(size(map,1)/length(pvalues)):size(map,1),:); % divide the colormap into number of needed levels.
newmap(1,:)=[0.5 0.5 0.5];
% colorbar tick labels
cbLables=cell(length(pvalues)+1,1);
cbLables{1}='0';
for pp=1:length(pvalues)
    cbLables{pp+1}=num2str(pvalues(pp));
end
colormap(newmap);
imagesc(sigShades_dataset1(:,:,length(pvalues)));
h=colorbar;
set(h,'YTickLabel',cbLables)
print('-dtiff ','-r72','dataset1sigShading')
close(gcf);

tfce_sig_dataset1= obstfcedataset1> v1(round((1-p)*nboot));
dataset1mapfinal=tfce_sig_dataset1.*obstmapdataset1;

if ~isempty(dataset2)
    all_tfce_sig_dataset2=zeros(xSize,ySize,length(pvalues));
    all_tfce_sig_contrast=zeros(xSize,ySize,length(pvalues));
    
    for ii=1:length(pvalues)
        thisp=pvalues(ii);
        all_tfce_sig_dataset2(:,:,ii)= obstfcedataset2> v2(round((1-thisp)*nboot));
        all_tfce_sig_contrast_pos= obstfcecontrastpos> vpos(round((1-thisp/2)*nboot));
        all_tfce_sig_contrast_neg= obstfcecontrastneg> vneg(round((1-thisp/2)*nboot));
        all_tfce_sig_contrast(:,:,ii)=(all_tfce_sig_contrast_pos+all_tfce_sig_contrast_neg)>0;
    end
    
    % shading significance map of dataset2
    sigShades_dataset2=zeros(xSize,ySize,length(pvalues));
    sigShades_dataset2(:,:,1)=all_tfce_sig_dataset2(:,:,1);
    for ii=2:length(pvalues)
        sigShades_dataset2(:,:,ii)=(ones(xSize,ySize)-all_tfce_sig_dataset2(:,:,ii-1)).*all_tfce_sig_dataset2(:,:,ii)*ii+sigShades_dataset2(:,:,ii-1);
    end
    
    figure; colormap(newmap);
    imagesc(sigShades_dataset2(:,:,length(pvalues)));
    h=colorbar;
    set(h,'YTickLabel',cbLables)
    print('-dtiff ','-r72','dataset2sigShading')
    close(gcf);
    
    
    % shading significance map of contrast map
    sigShades_contrast=zeros(xSize,ySize,length(pvalues));
    sigShades_contrast(:,:,1)=all_tfce_sig_contrast(:,:,1);
    for ii=2:length(pvalues)
        sigShades_contrast(:,:,ii)=(ones(xSize,ySize)-all_tfce_sig_contrast(:,:,ii-1)).*all_tfce_sig_contrast(:,:,ii)*ii+sigShades_contrast(:,:,ii-1);
    end
    
    figure; colormap(newmap);
    imagesc(sigShades_contrast(:,:,length(pvalues)));
    h=colorbar;
    set(h,'YTickLabel',cbLables)
    print('-dtiff ','-r72','contrastSigShading')
    close(gcf);
    
    % the final significance map for the given p level
    tfce_sig_dataset2= obstfcedataset2> v2(round((1-p)*nboot));
    tfce_sig_contrast_pos= obstfcecontrastpos> vpos(round((1-p/2)*nboot));
    tfce_sig_contrast_neg= obstfcecontrastneg> vneg(round((1-p/2)*nboot));
    tfce_sig_contrast=(tfce_sig_contrast_pos+tfce_sig_contrast_neg)>0;
    dataset2mapfinal=tfce_sig_dataset2.*obstmapdataset2;
    contrastmapfinal=tfce_sig_contrast.*obstmapcontrast;
    
end


%% effect size
% if ~isempty(dataset2)
%     % Effect sizes
%     if sum(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0)>0)>0
%         d_area1=(mean(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0))-mean(dataset2map(all_tfce_sig_contrast.*obstmapcontrast>0)))/sqrt((std(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0))^2 + std(dataset2map(all_tfce_sig_contrast.*obstmapcontrast>0))^2)/2);
%     end
%     if sum(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0)>0)>0
%         d_area2=(mean(dataset1map(all_tfce_sig_contrast.*obstmapcontrast<0))-mean(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0)))/sqrt((std(dataset1map(all_tfce_sig_contrast.*obstmapcontrast<0))^2 + std(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0))^2)/2);
%     end
%     if sum(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0)>0)==0 && sum(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0)>0)>0
%         Zscore=[mean(dataset1map(all_tfce_sig_contrast>0)) mean(dataset2map(all_tfce_sig_contrast>0))];
%         cohend=[d_area2];
%     elseif sum(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0)>0)==0 && sum(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0)>0)>0
%         Zscore=[mean(dataset1map(all_tfce_sig_contrast>0)) mean(dataset2map(all_tfce_sig_contrast>0))];
%         cohend=[d_area1];
%     elseif sum(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0)>0)==0 && sum(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0)>0)==0
%         Zscore=[];
%         cohend=[];
%     else
%         Zscore=[mean(dataset1map(all_tfce_sig_contrast.*obstmapcontrast>0)) mean(dataset1map(all_tfce_sig_contrast.*obstmapcontrast<0)) mean(dataset2map(all_tfce_sig_contrast.*obstmapcontrast>0)) mean(dataset2map(all_tfce_sig_contrast.*obstmapcontrast<0))];
%         cohend=[abs(d_area1) abs(d_area2)];
%     end
%     dlmwrite('Zscore.txt', Zscore, 'delimiter', '\t', 'precision', 7);
%     dlmwrite('cohend.txt', cohend, 'delimiter', '\t', 'precision', 7);
% end

%% Make figure with background (if specified) and significant areas
% Export the group and difference maps in tiff format
dataset1mapfinal(isnan(dataset1mapfinal))=0;
if isfield(cfg,'scaledownup')
    scaledown=cfg.scaledownup(1);
    scaleup=cfg.scaledownup(2);
    figure, imagesc(dataset1mapfinal,[scaledown scaleup]), colorbar
else
    scaledown=[];
    scaleup=[];
    figure, imagesc(dataset1mapfinal), colorbar
end
print('-dtiff ','-r72','dataset1statmap')
close(gcf);

if isempty(dataset2)==0
    dataset2mapfinal(isnan(dataset2mapfinal))=0;
    contrastmapfinal(isnan(contrastmapfinal))=0;
    if isfield(cfg,'scaledownup')
        scaledown=cfg.scaledownup(1);
        scaleup=cfg.scaledownup(2);
        figure, imagesc(dataset2mapfinal,[scaledown scaleup]), colorbar
    else
        scaledown=[];
        scaleup=[];
        figure, imagesc(dataset2mapfinal), colorbar
    end
    print('-dtiff ','-r72','dataset2statmap')
    close(gcf);
    if isfield(cfg,'scaledownup')
        scaledown=cfg.scaledownup(1);
        scaleup=cfg.scaledownup(2);
        figure, imagesc(contrastmapfinal,[scaledown scaleup]), colorbar
    else
        scaledown=[];
        scaleup=[];
        figure, imagesc(contrastmapfinal), colorbar
    end
    print('-dtiff ','-r72','contraststatmap')
    close(gcf);
end

dataset1cmat_temp=colormap;
dataset1pic=indtorgb(dataset1mapfinal,scaledown,scaleup,dataset1cmat_temp);

if isempty(dataset2)==0
    dataset2cmat_temp=colormap;
    dataset2pic=indtorgb(dataset2mapfinal,scaledown,scaleup,dataset2cmat_temp);
    Zdiffcmat_temp=colormap;
    if nansum(contrastmapfinal(:))==0
        contrastmapfinal(:,:)=0;
        diffpic(:,:,1)=contrastmapfinal;
        diffpic(:,:,2)=contrastmapfinal;
        diffpic(:,:,3)=contrastmapfinal;
    else
        diffpic=indtorgb(contrastmapfinal,scaledown,scaleup,Zdiffcmat_temp);
    end
end
close all;

% add background if a background picture is specified
if isfield(cfg,'backgroundfile')
    % open the background picture
    imbackground = double(imread(sprintf(backgroundfile)))/255;
    % add maps to background
    if ndims(imbackground)==2
        im3D = zeros(xSize, ySize, 3);
        im3D(:, :, 1) = imbackground;
        im3D(:, :, 2) = imbackground;
        im3D(:, :, 3) = imbackground;
    elseif ndims(imbackground)==3
        im3D = zeros(xSize, ySize, 3);
        im3D(:, :, 1) = imbackground(:, :, 1);
        im3D(:, :, 2) = imbackground(:, :, 2);
        im3D(:, :, 3) = imbackground(:, :, 3);
    end
    dataset1pic2 = im3D * transpim + dataset1pic * transpmap;
    for ii=1:xSize
        for jj=1:ySize
            if all_tfce_sig_dataset1(ii,jj)==0
                dataset1pic2(ii,jj,:)=im3D(ii,jj,:);
            end
        end
    end
    if isempty(dataset2)==0
        diffpic2 = im3D * transpim + diffpic * transpmap;
        for ii=1:xSize
            for jj=1:ySize
                if all_tfce_sig_contrast(ii,jj)==0
                    diffpic2(ii,jj,:)=im3D(ii,jj,:);
                end
            end
        end
        dataset2pic2 = im3D * transpim + dataset2pic * transpmap;
        for ii=1:xSize
            for jj=1:ySize
                if all_tfce_sig_dataset2(ii,jj)==0
                    dataset2pic2(ii,jj,:)=im3D(ii,jj,:);
                end
            end
        end
    end
else
    dataset1pic2 = dataset1pic;
    if isempty(dataset2)==0
        dataset2pic2 = dataset2pic;
        diffpic2 = diffpic;
    end
end

% save
dataset1picedge= dataset1pic2 .*255;
name = sprintf('dataset1picedge.tiff');
imwrite(uint8(dataset1picedge), name)
if isempty(dataset2)==0
    diffpicedge= diffpic2.*255;
    name = sprintf('diffpicedge.tiff');
    imwrite(uint8(diffpicedge), name)
    
    dataset2picedge= dataset2pic2.*255;
    name = sprintf('dataset2picedge.tiff');
    imwrite(uint8(dataset2picedge), name)
end

%% Eye-tracking measures in the significant areas and optional independent random split-data effect sizes for each of the areas significantly different between datasets (according to the difference matrix)
if isempty(dataset2)==0
    % number of significant pixels in single and difference matrices
    nbpixeltotal=xSize*ySize;
    poscluster=contrastmapfinal>0;
    negcluster=contrastmapfinal<0;
    nbsignifpixelarea1=sum(sum(poscluster));
    nbsignifpixelarea2=sum(sum(negcluster));
    nbpixelother=nbpixeltotal-nbsignifpixelarea1-nbsignifpixelarea2;
    for datafilenumber=1:length(datatotal)
        summary=[];
        durationarea1=[];
        durationarea2=[];
        meanfixdurarea1=[];
        meanfixdurarea2=[];
        meanfixdurrest=[];
        pathlengtharea1=[];
        pathlengtharea2=[];
        pathlenthrest=[];
        totfixdurarea1=[];
        totfixdurarea2=[];
        totfixdurrest=[];
        numfixarea1=[];
        numfixarea2=[];
        numfixrest=[];
        
        datatoload=[inputPath '\data' num2str(datatotal(datafilenumber))];
        load(datatoload); % The name of the matrix is 'summary'
        [nbfix nbvariables]=size(summary);
        
        nbitem=0;
        nbfixarea1=0;
        nbfixarea2=0;
        nbfixrest=0;
        cumuldurationarea1=0;
        cumuldurationarea2=0;
        cumuldurationrest=0;
        cumulsaccadearea1=0;
        cumulsaccadearea2=0;
        cumulsaccaderest=0;
        coordX = round(summary(1, columny));
        coordY = round(summary(1, columnx));
        if coordX<xSize && coordY<ySize && coordX>0 && coordY>0
            if poscluster(coordX(1),coordY(1))~=0
                cumuldurationarea1= cumuldurationarea1+ summary(1, columnduration);
                nbfixarea1=nbfixarea1+1;
            elseif negcluster(coordX(1),coordY(1))~=0
                cumuldurationarea2= cumuldurationarea2+ summary(1, columnduration);
                nbfixarea2=nbfixarea2+1;
            else
                cumuldurationrest=cumuldurationrest+summary(1, columnduration);
                nbfixrest=nbfixrest+1;
            end
        end
        for fix=2:nbfix
            coordX = round(summary(fix, columny));
            coordY = round(summary(fix, columnx));
            
            if coordX>0 && coordY>0 && coordX<xSize && coordY<ySize
                if summary(fix,columnitem)==summary(fix-1,columnitem)
                    if poscluster(coordX,coordY)~=0
                        cumuldurationarea1= cumuldurationarea1+ summary(fix, columnduration);
                        nbfixarea1=nbfixarea1+1;
                        cumulsaccadearea1=cumulsaccadearea1+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
                    elseif negcluster(coordX,coordY)~=0
                        cumuldurationarea2= cumuldurationarea2+ summary(fix, columnduration);
                        nbfixarea2=nbfixarea2+1;
                        cumulsaccadearea2=cumulsaccadearea2+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
                    else
                        cumuldurationrest=cumuldurationrest+summary(fix, columnduration);
                        nbfixrest=nbfixrest+1;
                        cumulsaccaderest=cumulsaccaderest+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
                    end
                elseif summary(fix,columnitem)~=summary(fix-1,columnitem)
                    nbitem=nbitem+1;
                    meanfixdurarea1(nbitem)=cumuldurationarea1/nbfixarea1;
                    meanfixdurarea2(nbitem)=cumuldurationarea2/nbfixarea2;
                    meanfixdurrest(nbitem)=cumuldurationrest/nbfixrest;
                    pathlengtharea1(nbitem)=cumulsaccadearea1;
                    pathlengtharea2(nbitem)=cumulsaccadearea2;
                    pathlenthrest(nbitem)=cumulsaccaderest;
                    totfixdurarea1(nbitem)=cumuldurationarea1;
                    totfixdurarea2(nbitem)=cumuldurationarea2;
                    totfixdurrest(nbitem)=cumuldurationrest;
                    numfixarea1(nbitem)=nbfixarea1;
                    numfixarea2(nbitem)=nbfixarea2;
                    numfixrest(nbitem)=nbfixrest;
                    
                    nbfixarea1=0;
                    nbfixarea2=0;
                    nbfixrest=0;
                    cumuldurationarea1=0;
                    cumuldurationarea2=0;
                    cumuldurationrest=0;
                    cumulsaccadearea1=0;
                    cumulsaccadearea2=0;
                    cumulsaccaderest=0;
                    coordX = round(summary(1, columny));
                    coordY = round(summary(1, columnx));
                    if coordX>0 && coordY>0 && coordX<xSize && coordY<ySize
                        if poscluster(coordX,coordY)~=0
                            cumuldurationarea1= cumuldurationarea1+ summary(fix, columnduration);
                            nbfixarea1=nbfixarea1+1;
                        elseif negcluster(coordX,coordY)~=0
                            cumuldurationarea2= cumuldurationarea2+ summary(fix, columnduration);
                            nbfixarea2=nbfixarea2+1;
                        else
                            cumuldurationrest=cumuldurationrest+summary(fix, columnduration);
                            nbfixrest=nbfixrest+1;
                        end
                    end
                end
            end
        end
        
        coordX = round(summary(fix, columny));
        coordY = round(summary(fix, columnx));
        if coordX>0 && coordY>0 && coordX<xSize && coordY<ySize
            if poscluster(coordX,coordY)~=0
                cumuldurationarea1= cumuldurationarea1+ summary(fix, columnduration);
                nbfixarea1=nbfixarea1+1;
                cumulsaccadearea1=cumulsaccadearea1+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
            elseif negcluster(coordX,coordY)~=0
                cumuldurationarea2= cumuldurationarea2+ summary(fix, columnduration);
                nbfixarea2=nbfixarea2+1;
                cumulsaccadearea2=cumulsaccadearea2+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
            else
                cumuldurationrest=cumuldurationrest+summary(fix, columnduration);
                nbfixrest=nbfixrest+1;
                cumulsaccaderest=cumulsaccaderest+sqrt((summary(fix,columnx)-summary(fix-1,columnx)).^2+(summary(fix,columny)-summary(fix-1,columny)).^2);
            end
        end
        
        nbitem=nbitem+1;
        meanfixdurarea1(nbitem)=cumuldurationarea1/nbfixarea1;
        meanfixdurarea2(nbitem)=cumuldurationarea2/nbfixarea2;
        meanfixdurrest(nbitem)=cumuldurationrest/nbfixrest;
        pathlengtharea1(nbitem)=cumulsaccadearea1;
        pathlengtharea2(nbitem)=cumulsaccadearea2;
        pathlengthrest(nbitem)=cumulsaccaderest;
        totfixdurarea1(nbitem)=cumuldurationarea1;
        totfixdurarea2(nbitem)=cumuldurationarea2;
        totfixdurrest(nbitem)=cumuldurationrest;
        numfixarea1(nbitem)=nbfixarea1;
        numfixarea2(nbitem)=nbfixarea2;
        numfixrest(nbitem)=nbfixrest;
        
        mediandurationarea1(datafilenumber)=nanmedian(meanfixdurarea1);
        mediandurationarea2(datafilenumber)=nanmedian(meanfixdurarea2);
        mediandurationrest(datafilenumber)=nanmedian(meanfixdurrest);
        totalduration(datafilenumber)=mediandurationarea1(datafilenumber)+mediandurationarea2(datafilenumber)+mediandurationrest(datafilenumber);
        relativedurationarea1(datafilenumber)=(mediandurationarea1(datafilenumber)/nbsignifpixelarea1)/(totalduration(datafilenumber)/nbpixeltotal);
        relativedurationarea2(datafilenumber)=(mediandurationarea2(datafilenumber)/nbsignifpixelarea2)/(totalduration(datafilenumber)/nbpixeltotal);
        relativedurationrest(datafilenumber)=(mediandurationrest(datafilenumber)/nbpixelother)/(totalduration(datafilenumber)/nbpixeltotal);
        
        fixdurarea1sbj(datafilenumber)=nanmean(meanfixdurarea1');
        fixdurarea2sbj(datafilenumber)=nanmean(meanfixdurarea2');
        fixdurrestsbj(datafilenumber)=nanmean(meanfixdurrest');
        pathlengtharea1sbj(datafilenumber)=nanmean(pathlengtharea1');
        pathlengtharea2sbj(datafilenumber)=nanmean(pathlengtharea2');
        pathlengthrestsbj(datafilenumber)=nanmean(pathlengthrest');
        totfixdurarea1sbj(datafilenumber)=nanmean(totfixdurarea1');
        totfixdurarea2sbj(datafilenumber)=nanmean(totfixdurarea2');
        totfixdurrestsbj(datafilenumber)=nanmean(totfixdurrest');
        numfixarea1sbj(datafilenumber)=nanmean(numfixarea1');
        numfixarea2sbj(datafilenumber)=nanmean(numfixarea2');
        numfixrestsbj(datafilenumber)=nanmean(numfixrest');
    end
    
    nbdatapointdataset1=0;
    nbdatapointdataset2=0;
    for ii=1:length(datatotal)
        if ismember(datatotal(ii),dataset1)
            nbdatapointdataset1=nbdatapointdataset1+1;
            if sum(poscluster(:))~=0
                fixdurarea1dataset1(nbdatapointdataset1)=fixdurarea1sbj(ii);
                pathlengtharea1dataset1(nbdatapointdataset1)=pathlengtharea1sbj(ii);
                totfixdurarea1dataset1(nbdatapointdataset1)=totfixdurarea1sbj(ii);
                numfixarea1dataset1(nbdatapointdataset1)=numfixarea1sbj(ii);
            end
            if sum(poscluster(:))~=0 || sum(negcluster(:))~=0
                totfixdurrestdataset1(nbdatapointdataset1)=totfixdurrestsbj(ii);
                numfixrestdataset1(nbdatapointdataset1)=numfixrestsbj(ii);
                fixdurrestdataset1(nbdatapointdataset1)=fixdurrestsbj(ii);
                pathlengthrestdataset1(nbdatapointdataset1)=pathlengthrestsbj(ii);
            end
            if sum(negcluster(:))~=0
                fixdurarea2dataset1(nbdatapointdataset1)=fixdurarea2sbj(ii);
                pathlengtharea2dataset1(nbdatapointdataset1)=pathlengtharea2sbj(ii);
                totfixdurarea2dataset1(nbdatapointdataset1)=totfixdurarea2sbj(ii);
                numfixarea2dataset1(nbdatapointdataset1)=numfixarea2sbj(ii);
            end
        elseif ismember(datatotal(ii),dataset2)
            nbdatapointdataset2=nbdatapointdataset2+1;
            if sum(poscluster(:))~=0
                fixdurarea1dataset2(nbdatapointdataset2)=fixdurarea1sbj(ii);
                pathlengtharea1dataset2(nbdatapointdataset2)=pathlengtharea1sbj(ii);
                totfixdurarea1dataset2(nbdatapointdataset2)=totfixdurarea1sbj(ii);
                numfixarea1dataset2(nbdatapointdataset2)=numfixarea1sbj(ii);
            end
            if sum(negcluster(:))~=0
                fixdurarea2dataset2(nbdatapointdataset2)=fixdurarea2sbj(ii);
                pathlengtharea2dataset2(nbdatapointdataset2)=pathlengtharea2sbj(ii);
                totfixdurarea2dataset2(nbdatapointdataset2)=totfixdurarea2sbj(ii);
                numfixarea2dataset2(nbdatapointdataset2)=numfixarea2sbj(ii);
            end
            if sum(poscluster(:))~=0 || sum(negcluster(:))~=0
                fixdurrestdataset2(nbdatapointdataset2)=fixdurrestsbj(ii);
                pathlengthrestdataset2(nbdatapointdataset2)=pathlengthrestsbj(ii);
                totfixdurrestdataset2(nbdatapointdataset2)=totfixdurrestsbj(ii);
                numfixrestdataset2(nbdatapointdataset2)=numfixrestsbj(ii);
            end
        end
    end
    
    if sum(poscluster(:))==0 && sum(negcluster(:))~=0
        eyeareadataset1=[fixdurarea2dataset1' fixdurrestdataset1' pathlengtharea2dataset1' pathlengthrestdataset1' totfixdurarea2dataset1' totfixdurrestdataset1' numfixarea2dataset1' numfixrestdataset1'];
        eyeareadataset2=[fixdurarea2dataset2' fixdurrestdataset2' pathlengtharea2dataset2' pathlengthrestdataset2' totfixdurarea2dataset2' totfixdurrestdataset2' numfixarea2dataset2' numfixrestdataset2'];
    elseif sum(negcluster(:))==0 && sum(poscluster(:))~=0
        eyeareadataset1=[fixdurarea1dataset1' fixdurrestdataset1' pathlengtharea1dataset1' pathlengthrestdataset1' totfixdurarea1dataset1' totfixdurrestdataset1' numfixarea1dataset1' numfixrestdataset1'];
        eyeareadataset2=[fixdurarea1dataset2' fixdurrestdataset2' pathlengtharea1dataset2' pathlengthrestdataset2' totfixdurarea1dataset2' totfixdurrestdataset2' numfixarea1dataset2' numfixrestdataset2'];
    elseif sum(poscluster(:))==0 && sum(negcluster(:))==0
        eyeareadataset1=[];
        eyeareadataset2=[];
    else
        eyeareadataset1=[fixdurarea1dataset1' fixdurarea2dataset1' fixdurrestdataset1' pathlengtharea1dataset1' pathlengtharea2dataset1' pathlengthrestdataset1' totfixdurarea1dataset1' totfixdurarea2dataset1' totfixdurrestdataset1' numfixarea1dataset1' numfixarea2dataset1' numfixrestdataset1'];
        eyeareadataset2=[fixdurarea1dataset2' fixdurarea2dataset2' fixdurrestdataset2' pathlengtharea1dataset2' pathlengtharea2dataset2' pathlengthrestdataset2' totfixdurarea1dataset2' totfixdurarea2dataset2' totfixdurrestdataset2' numfixarea1dataset2' numfixarea2dataset2' numfixrestdataset2'];
    end
    
    nametosave = [outputPath '\eyeareadataset1.txt'];
    dlmwrite(nametosave, eyeareadataset1, 'delimiter', '\t', 'precision', 7);
    nametosave = [outputPath '\eyeareadataset2.txt'];
    dlmwrite(nametosave, eyeareadataset2, 'delimiter', '\t', 'precision', 7);
end

fclose(imapreport);
