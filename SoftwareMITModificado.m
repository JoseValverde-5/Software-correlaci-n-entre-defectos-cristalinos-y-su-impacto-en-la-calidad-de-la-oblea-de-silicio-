% Dislocation counting code for scanner images
% Copyright (C) 2013  Massachusetts Institute of Technology

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


%
%--------------------------Directions-------------------------
%Scan sample. Save in format supported by imread, .png seems to work well.
%Open image in photo editor and crop image to ROI (recommend using a corner to help with 
%registration of microscope and scanner images, however absolute edge should be cropped out of image to avoid artifacts).
%Measure the area of representative single isolated dislocations to set limits on min and max dislocation sizes in [px^2].
%Enter file path to image for Im variable
%Enter appropriate values in maxSingleDislocationArea and minSingleDislocationArea
%Enter remaining inputs shown below

close all
clear
tic
fprintf('\n\n');
disp('Dislocation Counter-----------------------------------------------------');
warning off;

%--------------------------Inputs-------------------------
I = imread('Threshold.tif'); %Image file name

imageHeight=4;                                  %Length of sample in x-direction [cm]
imageWidth=4;                                   %Length of sample in y-direction [cm]
resolutionX=.75;                          		%Desired mapping resolution in x-direction [mm]
resolutionY=.75;                          		%Desired mapping resolution in y-direction [mm]
minSingleDislocationArea=5;          			%Min observed single dislocation size [px^2] Recommendations: 40 for 12800 DPI, 30 for 9600 DPI, 15 for 4800 DPI, 5 for 1200 DPI, 2 for 600 DPI
threshold=0;                        			%Threshold, input 0 to 1, -1 for auto. Corresponds in photo editor to x/[max color depth: e.g. 255 for 16 bits]
plotFlag=1;                          			%0 for no plots, 1 for yes
trimLength=0;                       			%[cm] Edge of sample often requires cropping to avoid artifacts. Can be done in image processing software or with this part of code
calibrationFlag=1;                      		%1 for calibration completed, so dislocation densities are displayed, 0 for no calibration, just pixel counts are displayed
calibrationSlope = .148;                   		%Slope from linear fit of micrograph/scanner image comparison
calibrationYint = -1;                   		%Y-intercept from linear fit of micrograph/scanner image comparison

%--------------------------Filters-------------------------
if plotFlag == 1
 figure('name','Orginal Image'), imshow(I)                         %Shows original image
end
bwSource = im2bw(I,threshold);          %treshold aunque lo tengamos a 0, porque en este programa se puede realizar                            %Thresholds image convierte la imagen en binario basandose en threshold
clear I;                                                            %Frees up memory
bw=1-bwSource;                                                      %Inverts Image
%Filters out small dots
bwFiltered = bwareaopen(bw,minSingleDislocationArea); %Filters out dots that are too small to be dislocations
bwFilteredImage=1-bwFiltered;                                       %Inverts filtered image for user viewing            
if plotFlag == 1
    figure('name','Filtered Image'), imshow(bwFilteredImage)        %Shows filtered Image in positive
    figure('name','Inverted Image'), imshow(bwFiltered)             %Shows Inverted Image
end
bw=bwFiltered;
clear bwFilteredImage;                                              %Frees up memory

%--------------------------Calculations-------------------------
%Gets Size
[sizeY,sizeX]=size(bwSource);                                       %Gets size of input image for later use

%Loops through analysis to give data for segments of image as selected by user
Xsteps=round(10/resolutionX*imageWidth);      						%Number of segments in x-direction, dependent on resolution for mapping
Ysteps=round(10/resolutionY*imageHeight);     						%Number of segments in y-direction dependent on resolution for mapping
fprintf('\n');
loops=0;
pixelCount = zeros(Xsteps,Ysteps);
 for xLoopCounter=1:Xsteps
 	for yLoopCounter=1:Ysteps
        loops=loops+1;

%Determines boundaries [in px] for cropping on each run
xmin=floor((sizeX-1)/Xsteps)*xLoopCounter-floor((sizeX-1)/Xsteps)+1; % This MATLAB function rounds each element of X to the nearest integer less than
    %or equal to that element.
ymin=floor((sizeY-1)/Ysteps)*yLoopCounter-floor((sizeY-1)/Ysteps)+1;
xmax=floor((sizeX-1)/Xsteps)*xLoopCounter;
ymax=floor((sizeY-1)/Ysteps)*yLoopCounter;

%Determine size of trimming area from edge region
cropLengthY=round(sizeY*trimLength*imageHeight);%Redondear al decimal o entero más cercano
cropLengthX=round(sizeX*trimLength*imageWidth);
cornerTrim=(xmax-xmin)*cropLengthY+cropLengthX*(ymax-ymin-cropLengthY);
edgeTrimTop=(xmax-xmin)*cropLengthY;
edgeTrimSide=(ymax-ymin)*cropLengthX;

%Isolates individual segment, trims edges by amount specified by user, calculates total number of dislocated pixels
if xLoopCounter==1 && yLoopCounter==1
    bwSegmented=bwFiltered([ymin+cropLengthY:ymax],[xmin+cropLengthX:xmax]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(cornerTrim+areaSegmented)/areaSegmented;
elseif xLoopCounter==1 && yLoopCounter==Ysteps
    bwSegmented=bwFiltered([ymin:ymax-cropLengthY],[xmin+cropLengthX:xmax]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(cornerTrim+areaSegmented)/areaSegmented;
elseif xLoopCounter==Xsteps && yLoopCounter==1
    bwSegmented=bwFiltered([ymin+cropLengthY:ymax],[xmin:xmax-cropLengthX]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(cornerTrim+areaSegmented)/areaSegmented;
elseif xLoopCounter==Xsteps && yLoopCounter==Ysteps
    bwSegmented=bwFiltered([ymin:ymax-cropLengthY],[xmin:xmax-cropLengthX]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(cornerTrim+areaSegmented)/areaSegmented;
elseif xLoopCounter==1
    bwSegmented=bwFiltered([ymin:ymax],[xmin+cropLengthX:xmax]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(edgeTrimSide+areaSegmented)/areaSegmented;
elseif yLoopCounter==1
    bwSegmented=bwFiltered([ymin+cropLengthY:ymax],[xmin:xmax]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(edgeTrimTop+areaSegmented)/areaSegmented;
elseif yLoopCounter==Ysteps
    bwSegmented=bwFiltered([ymin:ymax-cropLengthY],[xmin+cropLengthX:xmax]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(edgeTrimTop+areaSegmented)/areaSegmented;
elseif xLoopCounter==Xsteps
    bwSegmented=bwFiltered([ymin:ymax],[xmin:xmax-cropLengthX]);
    areaSegmented=size(bwSegmented,1)*size(bwSegmented,2);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented))*(edgeTrimSide+areaSegmented)/areaSegmented;
else 
    bwSegmented=bwFiltered([ymin:ymax],[xmin:xmax]);
    pixelCount(yLoopCounter,xLoopCounter) = sum(sum(bwSegmented));
end
	end
 end	

%--------------------------Outputs-------------------------
%Displays matrix of counted etch pit pixels
disp('Total number of etch pit pixels');
pixelCount

if calibrationFlag==1
    %Displays matrix of dislocation densities
    dislocationDensity =(pixelCount+calibrationYint)/calibrationSlope/(imageWidth/Xsteps*imageHeight/Ysteps)
   
    beep;
    

    %Creates dislocation density map
    densityMap=dislocationDensity;
    filter=4e6>densityMap&1e4<densityMap;   						%Creates matrix to pick out dislocation densities within the detection limits of measurement
    setMax=4e6<=densityMap;     									%Creates matrix to filter out dislocation densities greater than the upper detection limit
    setMin=1e4>=densityMap;    		 								%Creates matrix to filter out dislocation densities less than lower detection limit
    densityMap=filter*densityMap+(4e6*setMax)+(1e4*setMin);    	%Sets all dislocation densities above upper limit to upper limit and all below lower limit to lower limit
    rescale=log10(densityMap);  									% muestra los datos de la matriz C como una imagen que utiliza el rango completo de colores en el colores. Cada elemento de C especifica el color de 1 píxel de la imagen. Puts dislocation densities on a log scale
    grays=transpose([1:-1/128:0;1:-1/128:0;1:-1/128:0]);    		%Creates grayscale for dislocation density map
    %Creates contour plot of dislocation density vs position
    densityMap=dislocationDensity;
    filter=4e6>densityMap&1e3<densityMap;                           %Creates matrix to pick out dislocation densities within the detection limits of measurement
    setMax=4e6<=densityMap;                                         %Creates matrix to filter out dislocation densities greater than the upper detection limit
    setMin=1e3>=densityMap;                                         %Creates matrix to filter out dislocation densities less than lower detection limit
    densityMap=filter.*densityMap+(4e6*setMax)+(1e3*setMin);        %Sets all dislocation densities above upper limit to upper limit and all below lower limit to lower limit
    rescale=log10(densityMap);                                      % muestra los datos de la matriz C como una imagen que utiliza el rango completo de colores en el colores. Cada elemento de C especifica el color de 1 píxel de la imagen. Puts dislocation densities on a log scale
    grays=transpose([1:-1/128:0;1:-1/128:0;1:-1/128:0]);            %Creates grayscale for dislocation density map
    %Creates contour plot of dislocation density vs position
    figure('name','Dislocation Density Map'), imagesc(rescale)
    hold on   
    axis off
    axis image
    colormap(grays)
    colorbar('location','eastoutside','YTick',...
[min(min(rescale)) (3*min(min(rescale))+max(max(rescale)))/4 (min(min(rescale))+max(max(rescale)))/2 (min(min(rescale))+3*(max(max(rescale))))/4 max(max(rescale))],...
'YTickLabel',...
{num2str(10^min(min(rescale)),'%0.2G'),num2str(10^((3*min(min(rescale))+max(max(rescale)))/4),'%0.2G'),num2str(10^((min(min(rescale))+max(max(rescale)))/2),'%0.2G'),num2str(10^((min(min(rescale))+3*(max(max(rescale))))/4),'%0.2G'),num2str(10^max(max(rescale)),'%0.2G')})
    hold off
 
end
colormap(jet)% Gama de colores"jet"
saveas(gcf,'D1.png')
fprintf('Execution Time = %4.3f [sec]\n',toc);
beep;
