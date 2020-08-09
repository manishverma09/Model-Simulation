function [blueSoilPixels, redSoilPixels, nirSoilPixels, mirSoilPixels, ndviSoilPixels] = findPureSoilPixels (plotOrNot)
% THis function tries to locate pure soil pixels using MOD13 in the area around Sturt
% Plains site.

% output: Outputs blue, red, nir, and mir band of identified pixels

% arguments: whethere you want the function to produce plots or not. Useful
%           to get the plots out when you are examining but not once you are
%           satisfied with the result.

%           plotOrNot: 'Y' (for yes) or 'N' (for no);



%% I tried different methods to try and locate pure soil pixels. Criteria

%1. Red and inr scatter near soil line
%2. Absence of vegetation seasonality
%3. Low peak NDVI
%4. Inverse realtion between moisture and reflectance
%5. mir>nir>red>blue - red should be significantly higher than blue

% Finally, I am using the idea of locating pixels with smallest peak NDVI.
% To confirm that the identified pixels are likely to be pure-soil, I look
% at the seasonality of NDVI (it should not show vegetation growth),
% compare red and nir reflectance (the scatter should be close to 1:1line
% and should not show the distinct triangular pattern shown by vegetated
% pixels), verify that change in reflectance is largely because of
% moisture, and confirm that nir>red>blue - red should be significantly higher than blue.

%% Cd to data directory
cd /mod13/fiftyKm/GTiff

%% ########################
% we want to first use VI and find soil pixels
%###########################

%% Read VI data
fileList = dir('*NDVI*');     % dir('*EVI*')

trialData = geotiffread(fileList(1).name);
eviData = nan([size(trialData) length(fileList)]);

for i=1:numel(fileList)
    
    temp = geotiffread(fileList(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    
    eviData(:,:,i) = temp/10000;
end

%% Get date and year
fileList = dir('*NDVI*');

yearDoy = nan([length(fileList) 2]);

for i=1:numel(fileList)
    
    temp = strsplit(fileList(i).name, '.');
    
    temp = temp{2};
    yearVal = str2double(temp(2:5));
    doyVal = str2double(temp(6:8));
    
    yearDoy(i,:) = [yearVal, doyVal];
end

%% Get for the south hemisphere growing season

evi = eviData(:,:,15:end);

% find seasonal amplitude for each pixel and select the ones with smallest
% seasonal amplitude

% minEvi = min(evi, [], 3);
% maxEvi = max(evi, [], 3);
%
% eviAmp = maxEvi - minEvi;

%[row, col, index] = find(eviAmp<0.15);

%OR

% find the max of each pixel and select the ones with small maximum NDVI

eviSort = sort(evi, 3);   % ascending sort
maxEvi = eviSort(:,:,end-1:end);  % take the top three
maxEvi = nanmean(maxEvi, 3);

[row, col, ~] = find(maxEvi<0.32);

%% Plot for some and see

if strcmp(plotOrNot, 'Y') ==1
    
    figure;
    for i = 1:numel(row)
        
        plot(squeeze(evi(row(i), col(i), :)), '*-');
        pause;
        close all;
    end
end;
%% ###################################
%  Check via NIR and Red reflectance
%##################################

%% Read reflectance data
fileList1 = dir('*NIR*');
fileList2 = dir('*red*');
fileList3 = dir('*blue*');
fileList4 = dir('*MIR*');
fileList5 = dir('*composite_day_of_the_year*');

trialData = geotiffread(fileList1(1).name);

nirData = nan([size(trialData) length(fileList1)]);
redData = nan([size(trialData) length(fileList1)]);
doyData = nan([size(trialData) length(fileList1)]);
blueData = nan([size(trialData) length(fileList1)]);
mirData = nan([size(trialData) length(fileList1)]);


for i=1:numel(fileList1)
    
    temp = geotiffread(fileList1(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    nirData(:,:,i) = temp/10000;
    
    temp = geotiffread(fileList2(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    redData(:,:,i) = temp/10000;
    
    temp = geotiffread(fileList3(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    blueData(:,:,i) = temp/10000;
    
    temp = geotiffread(fileList4(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    mirData(:,:,i) = temp/10000;
    
    temp = geotiffread(fileList5(i).name);
    temp = double(temp);
    temp (temp <0) = nan;
    doyData(:,:,i) = temp;
end

%% Plot NIR and Red for the row and col found above using NDVI or EVI

if strcmp(plotOrNot, 'Y') ==1
    
    nirData1 = nirData(:,:,15:end);
    redData1 = redData(:,:,15:end);
    
    
    for i=1:numel(row)
        nir1 = squeeze(nirData1(row(i),col(i),:));
        red1 = squeeze(redData1(row(i),col(i),:));
        
        scatter(red1, nir1);
        box on;
        xlim([0.09 0.32])
        ylim([0.09 0.32])
        hold on;
        plot(xlim, ylim, 'k-')
        %hold all;
        pause;
        close all;
    end
    
    box on;
    xlim([0.09 0.32])
    ylim([0.09 0.32])
    hold on;
    plot(xlim, ylim, 'k-')
end

%%  Use diff between nir and red to find pure soil pixel
%
% nirMinusRed = nanmean(abs(nirData(:,:,15:end) - redData(:,:,15:end)), 3);
%
% nirData1 = nirData (:,:,15:end);
% redData1 = redData (:,:,15:end);
%
% [row, col, index1] = find(nirMinusRed<0.08);
%
% % plot nir and red and NDVI and check
%
% for i=1:numel(row)
%
%     ndviVal = (squeeze(nirData1(row(i), col(i), :)) - squeeze(redData1(row(i), col(i), :)));
%     ndviVal = ndviVal./(squeeze(nirData1(row(i), col(i), :)) + squeeze(redData1(row(i), col(i), :)));
%
%     subplot(1,2,1)
%     plot(ndviVal, '*-');
%
%     subplot(1,2,2)
%     scatter(squeeze(redData1(row(i), col(i), :)), squeeze(nirData1(row(i), col(i), :)))
%     box on;
%     xlim([0.06 0.26])
%     ylim([0.06 0.26])
%     hold on;
%     plot(xlim, ylim, 'k-')
%     pause;
%     close all;
% end

%% Plot and look for soil line
% for i=1:10:200
%     for j=1:10:300
%
%         xData = squeeze(redData(i,j,:));
%         yData = squeeze(nirData(i,j,:));
%
%         if sum(isnan(xData))~=numel(xData)
%
%
%         scatter(squeeze(redData(i,j,:)), squeeze(nirData(i,j,:)));
%
%         val1 = [min(squeeze(redData(i,j,:))), min(squeeze(nirData(i,j,:)))];
%         val1 = min(val1);
%
%         val2 = [max(squeeze(redData(i,j,:))), max(squeeze(nirData(i,j,:)))];
%         val2 = max(val2);
%
%         hold on;
%         xlim([val1 val2]);
%         ylim([val1 val2]);
%         plot(xlim, ylim, 'k-');
%
%         pause;
%         close all;
%         end
%     end
% end

%% Compare soil mositure and reflectance -

loadFileName = 'metAndRsData/SturtPlains_v12a_swInFparEvi';
load(loadFileName);

loadFileName = 'gppData/SturtPlains_v12a';
load(loadFileName);

index1= ismember(yearDoyTimeGpp(:,1), [2014, 2015]);
yearDoyTime = yearDoyTimeGpp(index1,1:3);

sm3cm = sm3cm(index1);

index1 = (abs(yearDoyTime(:,3)-0.6042)<0.0001);
yearDoyTime = yearDoyTime(index1,1:3);
sm3cm = sm3cm(index1);

yearDoyTime (yearDoyTime(:,1)==2015, 2) = yearDoyTime (yearDoyTime(:,1)==2015, 2) + 365;


% Get doy for modis reflectance of identified pixels
nirData = nirData(:,:,15:end);
redData = redData(:,:,15:end);
blueData = blueData(:,:,15:end);
mirData = mirData(:,:,15:end);
doyData = doyData(:,:,15:end);  % this is the day observations were recorded
yearDoy = yearDoy(15:end,:);    % the doy here is probably middle of the composite period

% Make days unique
index1 = yearDoy(:,1) ==2015;

if strcmp(plotOrNot, 'Y') ==1
    
    for i=1:numel(row)
        
        doyData1 = squeeze(doyData(row(i), col(i), :));
        doyData1(index1,1) = doyData1(index1,1)+365;
        
        
        index = ismember(yearDoyTime(:,2), doyData1);
        
        indexRev = ismember(doyData1, yearDoyTime(:,2));
        
        if isempty(indexRev) ==0
            
            sm3Val = sm3cm(index);
            nir1 = squeeze(nirData(row(i),col(i),:));
            red1 = squeeze(redData(row(i),col(i),:));
            
            subplot(1,2,1)
            scatter(sm3Val, nir1(indexRev))
            subplot(1,2,2)
            scatter(sm3Val, red1(indexRev));
            
            pause;
            close all;
        end
    end
end
%% Plot the four reflectance bands for the selected soil pixels and choose the ones that look the best

evi = eviData(:,:,15:end);

if strcmp(plotOrNot, 'Y') ==1
    for i=1:numel(row)
        
        
        nirVal = squeeze(nirData(row(i),col(i),:));
        redVal = squeeze(redData(row(i),col(i),:));
        blueVal = squeeze(blueData(row(i),col(i),:));
        mirVal = squeeze(mirData(row(i),col(i),:));
        eviVal = squeeze(evi(row(i),col(i),:));
        
        plot(nirVal, '*-');
        hold all;
        plot(redVal, '*-');
        plot(blueVal, '*-');
        plot(mirVal, '*-');
        plot(eviVal, 'X-');
        
        legend('nir', 'red', 'blue', 'mir', 'ndvi');
        pause;
        close all;
    end
end
%% Output the values

nirSoilPixels = nan(numel(row), size(nirData, 3));
redSoilPixels = nan(numel(row), size(nirData, 3));
blueSoilPixels = nan(numel(row), size(nirData, 3));
mirSoilPixels = nan(numel(row), size(nirData, 3));
ndviSoilPixels = nan(numel(row), size(nirData, 3));

for i=1:numel(row)
    
    
    nirSoilPixels(i,:) = squeeze(nirData(row(i),col(i),:));
    redSoilPixels(i,:) = squeeze(redData(row(i),col(i),:));
    blueSoilPixels(i,:) = squeeze(blueData(row(i),col(i),:));
    mirSoilPixels(i,:) = squeeze(mirData(row(i),col(i),:));
    ndviSoilPixels(i,:) = squeeze(evi(row(i),col(i),:));
    
end

