%% Threshold an image based on the gamma distribution fitting its histogram
% |Copyright 2017, Luca Della Santina|
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% This software is released under the terms of the GPL v3 software license
%
% |*GammaFitting fits an image stack with gamma function and masks the
% image to include only those pixels whose value exceeds a given a mount
% of times the standard deviation of the gamma fitting function.*|
%
% This program is especially useful to segment punctate labeling in
% immunofluorescence confocal images (i.e. synaptic proteins).
%
% Follows the procedure published in Note5/Fig.10 of the book chapter:
%
%  Using Fluorescent Markers to Estimate Synaptic Connectivity In Situ.
%  Hoon M, Sinha R, Okawa H. Methods Mol Biol. 2017;1538:293-320.
%
% <<Theory.png>>
%
% In addition, users choose whether a single SD threshold is
% selected and applied to the entire image stack or individual
% thresholds will be calculated for each Z-plane of the image stack.
% The latter allows to compensate for changes in the noise distribution
% due to uneven antibody penetration in the tissue or laser power correction
% applied during image acquisition.
% The distribution of applied thresholds is plotted as a function
% of Z-position in the image stack. The raw fitted values are smoothed (red
% line) to allow for errors of the nonlinear fitting estimators.
%
% * If only one threshold is applied to the entire image, the output will
%   be a distribution of the pixel intensities, with overlayed the gamma
%   fittinng curved used to calculate the threshold value
%
% <<SingleThreshold.png>>
%
% * If multiple thresholds are calculated, the output will be the
%   distribution of threshold values as a function of depth in the volume
%
% <<ThresholdChangeZaxis.png>>
%
% *Input:*
%
% * Source signam image (single-channel TIF image stack)
%
% *Output:*
%
% * Binary mask image (single-channel TIF image stack)
%
% *Dependencies:*
%
% * textprogressbar.m
%

% Load original image file
[FileName, PathName] = uigetfile('*.tif', 'Load a single channel tif file');
tmpImInfo = imfinfo([PathName FileName]);
tmpImSize = [tmpImInfo(1).Height tmpImInfo(1).Width length(tmpImInfo)];
textprogressbar(['Loading image stack from "' FileName '" ...']);
%disp(['Loading image stack from "' FileName '" ...']);
for i = 1:tmpImSize(3)
    tmpStack(:,:,i)=imread([PathName FileName], i);
    textprogressbar(100*i/tmpImSize(3)); % update progress bar
end
textprogressbar('DONE');

% Ask user for Options

tmpPrompt = {'How many SDs away from mode you want to set threshold? ',...
             'Fit entire stack with a single threshold? 0=no, 1=yes',...
             'Show debug info? 0=no, 1=warnings, 2=individual histograms',...
             'gamma-fit initial rise speed coefficient',...
             'gamma-fit initial decay speed coefficient'};
tmpAns = inputdlg(tmpPrompt, 'Threshold value',[1 40],{'4','0','0','0.1','50'});

tmpNumSDs      = str2num(tmpAns{1});   % number of gamma-fit standard devs.
fitEntireImage = str2num(tmpAns{2});   % fit mode (z-depth independency)
debug          = str2num(tmpAns{3});   % debug code active?
fitBeta0rise   = str2num(tmpAns{4});   % gamma-fit initial beta values
fitBeta0decay  = str2num(tmpAns{5});   % gamma-fit initial beta values

% Fit the image histogram with gamma function

if fitEntireImage == 1
    textprogressbar('Fitting the entire stack with a single gamma function: ');
else
    textprogressbar('Fitting each image along the z-axis with a different gamma function: ');
end
tmpThresh = zeros(1,tmpImSize(3));
tmpFitErrors = '';

for i = 1:tmpImSize(3)
    if debug==1 disp(i); end
    tmpI = tmpStack(:,:,i);              % extract the frame of interest
    IVectorized = reshape(tmpI, [1, size(tmpI,1)*size(tmpI,2)]);

    if fitEntireImage == 1
        IVectorized = reshape(tmpStack, [1, size(tmpStack,1)*size(tmpStack,2)*size(tmpStack,3)]);
    end

    if max(max(max(tmpStack))) < 256
        Edges = 0:1:max(max(max(tmpStack))); % 8-bit input image (max=255)
    else
        Edges = 0:44:max(max(max(tmpStack))); % for 12-bit or 16-bit images
    end
    [n, bin] = histc(IVectorized, Edges);
    n=smooth(n, 3, 'rloess');  % smooth histogram to get rid of local peaks
    n=n.';  % transpose the smoothed histogram back to a column vector


    if debug==2
        figure(); set(gcf, 'Color', [0 0 0]); set(gca, 'Color', [0 0 0]);
        hold on;
        stairs(Edges, n, 'w');
        set(gca, 'XColor', [1 1 1]); set(gca, 'YColor', [1 1 1]);
        xlabel('Intensity', 'Color', 'w');
        ylabel('Number of voxels', 'Color', 'w');
    end

    FitStartX = find(n==max(n))-1; % Find position of histogram max
    if max(max(max(tmpStack))) < 256
        FitStartX = FitStartX * 1;    % convert back to intensity values
    else
        FitStartX = FitStartX * 44;    % convert back to intensity values (LDS: fix this arbitrary 44 value)
    end
    FitStartX = FitStartX(end);    % If multiple maximum take the last

    Edges1=Edges(1:end-1); % Smooth the histrogram by averaging neighbors
    Edges2=Edges(2:end);
    EdgesMid = mean([Edges1;Edges2],1); % Edges is the x on the bin posizion on histogram
    n2 = n(1:end-1); % n2 is the histogram count of voxel at the "Edge" bin position
    FitStartXInd = find(EdgesMid>FitStartX, 1, 'first');
    EdgesMidFit = EdgesMid(FitStartXInd:end);
    n2Fit = n2(FitStartXInd:end);

    try
        %Fit gamma function with a modification that allows x shift.
            %beta0(1)=gamma distribution rising speed,
            %beta0(2)=gamma distribution falling speed,
        %beta0 = [max(n) fitBeta0rise fitBeta0decay FitStartX];
        beta0 = [fitBeta0rise fitBeta0decay];
        options = optimset('MaxIter',5000);
        if debug==0 warning('off','all'); end % Disable warnings from nlinfit when not debugging
        modelfun = @(beta,x)(max(n).*(x-FitStartX).^beta(1).*exp(-(x-FitStartX)./beta(2)));
        fitcoef = nlinfit(EdgesMidFit, n2Fit, modelfun, beta0, options);
        if debug==0 warning('on','all'); end %
        fitcoef = real(fitcoef); % taking real value works if complex number is returned
        fit = modelfun(fitcoef, EdgesMid);
        if debug==2 hold on, plot(EdgesMid(FitStartXInd:end), fit(FitStartXInd:end), 'r');end
        GammaFitMean(i) = (fitcoef(1)+1)*fitcoef(2)+FitStartX;
        GammaFitMode(i) = fitcoef(1)*fitcoef(2)+FitStartX;
        GammaFitVar(i) = (fitcoef(1)+1)*fitcoef(2)^2;
        GammaFitSD(i) = sqrt(abs(GammaFitVar(i))); %when fitting problems GammaFitVar is a negative number, try to fix better than taking absolute value
        tmpThresh(i) = GammaFitMode(i) + GammaFitSD(i) * tmpNumSDs(1);
    catch ME
        tmpFitErrors = [tmpFitErrors 'Error z-plane #' num2str(i) ', because:' ME.message '\n'];
    end
    textprogressbar(100*i/tmpImSize(3)); % update progress bar
    if fitEntireImage == 1
        break
    end
end
textprogressbar('DONE');
disp(sprintf(tmpFitErrors)); % Print errors happened during fitting on screen

% Plot fitting quality against z-axis or the entire stack histogram

if fitEntireImage == 0
    % Plot distribution of thresholds across z-axis to check fitting quality
    figure, plot(tmpThresh); hold on;
    title('blue=original fit, red=smoothed fit)');
    tmpThresh = smooth(tmpThresh, 10, 'rloess'); % Smooth thresholds to avoid spkes
    plot(tmpThresh, 'r');
    ylim([0 max(max(max(tmpStack)))]);
    ylabel('Threshold value'); xlabel('Image # along Z-axis');
else
    % Plot gamma-fit against the entire image stacks's histogram
    figure, stairs(Edges, n, 'k'), hold on;
    title('black=histogram, red=gamma fit, blue=threshold');
    plot([tmpThresh(1) tmpThresh(1)], [0 max(n)], 'b');
    plot(EdgesMid(FitStartXInd:end), fit(FitStartXInd:end), 'r');
    ylabel('Number of voxels'), xlabel('Intensity value');
    xlim([0 tmpThresh(1)*1.2]);
end

% Mask the original image stack with the gamma-fit thresholds
for i=1:numel(tmpNumSDs)
    % Recalculate tmpThresh using the selected number of SDs
    for j=1:tmpImSize(3)
        tmpThresh(j) = GammaFitMode(j) + GammaFitSD(j) * tmpNumSDs(i);
        if fitEntireImage == 1
            break
        end
    end
    tmpThresh = smooth(tmpThresh, 10, 'rloess'); % Smooth thresholds to avoid spkes

    textprogressbar(['Thresholding original image stack (' num2str(tmpNumSDs(i)) ' SD): ']);
    tmpStackMasked = uint8(tmpStack);
    for k = 1:tmpImSize(3)
        tmpI = tmpStack(:,:,k);                 % Get the frame of interest
        if fitEntireImage == 0
            tmpMask = tmpI>tmpThresh(k);        % Create mask from threshold
        else
            tmpMask = tmpI>tmpThresh(1);        % Create mask from threshold
        end
        tmpStackMasked(:,:,k)= uint8(tmpMask);  % Populate the 8-bit mask stack
        textprogressbar(100*k/tmpImSize(3));    % update progress bar
    end
    textprogressbar('DONE');

    % Save resulting mask stack into a TIF file (LZW compressed to save space)

    [FileName,PathName,FilterIndex] = uiputfile('*.tif','Name&Save the tif file.');
    textprogressbar(['Saving Mask into "' FileName '": ']);
    imwrite(tmpStackMasked(:,:,1), [PathName FileName], 'tif', 'compression', 'lzw'); %write first z plane
    if tmpImSize(3) > 1; %write the rest of the z planes
        for k=2:tmpImSize(3)
            imwrite(tmpStackMasked(:,:,k), [PathName FileName], 'tif', 'compression', 'lzw', 'WriteMode', 'append');
            textprogressbar(100*k/tmpImSize(3));    % update progress bar
        end
    end
    textprogressbar('DONE');
end
% Final cleanings
clear tmp* i FileName PathName fit* Fit* beta0 bin Edges* Filter* Gamma* model* n* options IVectorized debug;

%% Change log
% _*Version 2.5.1*             created on 2017-11-21 by Luca Della Santina_
%
%  + Reformatted documentation using proper MATLAB markup format
%  + Added Copyright info
%
% _*Version 2.5*               created on 2017-07-21 by Luca Della Santina_
%
%  + Allows inserting multiple SD threshold values to test which is best
%  % Output plot in single threshold shows properly scaled x-axis
%
% _*Version 2.4*               created on 2017-05-04 by Luca Della Santina_
%
%  % Image histogram distribution is smoothed (kern=3) before gamma fit
%     to prevent problem with images containing a lot of black voxels
%  % Single-gamma mode: Threshold is plotted (blue line) on the final plot
%
% _*Version 2.3*               created on 2017-04-08 by Luca Della Santina_
%
%  % Changed smoothing of the threshold curve from default moving averate
%    To the robust 2nd order polinomial rloess to remove outliers
%    happening when an error occurs during fitting
%  % Fixed typo in debug output text
%
% _*Version 2.2*               created on 2017-04-06 by Luca Della Santina_
%
%  + Errors during fitting loop are reported at the end with z-plane #
%  % Removed gamma function amplitude and start x from nlinfit beta0
%    as these parameters can be precisely calculated from the histogram
%    so they should not be further optimized (and changed) by fitting
%  % Histogram for 8-bit images is now computed and fitted at full
%    resolution to minimize fitting errors due to undersampling
%
% _*Version 2.1*               created on 2017-04-05 by Luca Della Santina_
%
%  + 8-bit input images are now recognized
%  + User is allowed to input custom initial paramenters for the gamma fit
%  + Progress bars displayed during each computationally-intensive step
%  % Output threshold vs z-depth plot now scales from 0 to image max
%  % Changed default initial gamma fit decay to 0.1 (more consistency)
%  % Changed gamma fit SD computation to not return error when poor fit
%
% _*Version 2.0*               created on 2017-04-03 by Luca Della Santina_
%
%  + User can choose whether to fit the entire stack's histogram with a
%    single gamma function or to fit each image along z-axis with different
%    gamma functions. The latter allows more effective masking of labelings
%    whose noise distribution changes along the z-axis.
%    These differences across z-axis usually arise due to:
%    . Penetration issues of the primary antibody
%    . Labeling more intense in some cell types than others (CtBP2)
%    . Laser correction along z-axis applied during image acquisition
%  + Thresholds are smoothed along z-axis to compensate errors of nlnfit
%    in some planes particularly hard to fit.
%  + Auto-detection of histogram peak (before value was user-requested)
%  + Threshold is asked to user via graphic dialog, proposing 4 as default
%  + Threshold distribution is plocltted against image number along z-axis
%  % Default value for beta0(1) = 1 instead of 0.5. (more consistent fits)
%  % FitStartX now matches histogram peak, before it fell 1 bin after peak
%
% _*Version 1.0*                                 created by Haruhisa Okawa_
%
%  + Fitting parameters are calculated once over the entire image stack
%  + User is asked to choose the intensity value of histogram's peak
%  + Output mask image is saved as uncompressed TIF stack
