function [meanOptimizedParameters, stdOptimizedParameters] =  compareSimulatedWithModis (N_value_dir, noOfSim)
% Compare PROSAIL results with MODIS observations and returns parameters
% that result in the minimum discrepency between simulated and MODIS
% observations

% output: 
%         meanOptimizedParameters: is the time series of PROSAIL parameters 
%         stdOptimizedParameters: standard deviation of PROSAIL parameters

% argument:
%         N_value_dir: is the directory where simulation results are saved.
%                      N value refers to the values of leaf structure parameter for PROSAIL
%                      simulation. Directories are named like
%                      'nOntPointFour', 'nOnePointFive' etc
%                      I have three directories nOnePointFour, nOnePointFive, and nOnePointSix 
%         noOfSim: This is based on deviation between hundreds of PROSAIL simulations and MODIS observations.
%                  It is the no of different prosail simulations you want to consider where the difference
%                  between simulated and MODIS reflectance is smaller than the rest of the simulations. Typical values is 100. 

%% Load MODIS data

load ('/workspace/reflectanceLandBand.mat');

%% Extract for 2014 AND 2015
%% First calculate mean and std of EVI - note that in the paper we are using only EVI

index = yearDoy(:,1) ==2014 | yearDoy(:,1) == 2015;
yearDoy = yearDoy(index,:);

band1 = reflectanceLandBand.Band1(6:8,6:8,index);
band2 = reflectanceLandBand.Band2(6:8,6:8,index);
band3 = reflectanceLandBand.Band3(6:8,6:8,index);

evi = (2.5)* (band2-band1)./ (band2 + 6*band1 - 7*band3 + 1);
ndvi = (band2-band1)./(band2+band1);

% 31 is day 241 in 2014 and 76 is day 233 in 2015

evi = (evi(:,:,31:76));    
ndvi = ndvi(:,:,31:76);
yearDoy = yearDoy(31:76,:);

smoothEvi = nan(3,3,46);
smoothNdvi = nan(3,3,46);
for i=1:3
    for j=1:3
    smoothEvi(i,j,:) = smooth(evi(i,j,:),3);
    smoothNdvi(i,j,:) = smooth(ndvi(i,j,:),3);

    end
end
    
meanEvi = nanmean(squeeze(nanmean(smoothEvi, 1)), 1);    
stdEvi = nanstd(squeeze(nanmean(smoothEvi, 1)), 1);   

meanNdvi = nanmean(squeeze(nanmean(smoothNdvi, 1)), 1);    
stdNdvi = nanstd(squeeze(nanmean(smoothNdvi, 1)), 1);     
    
%% Now get separate bands - you may have to chage the foolowing if you want std values
band1 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band1(6:8,6:8,index), 1)), 1));
band2 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band2(6:8,6:8,index), 1)), 1));
band3 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band3(6:8,6:8,index), 1)), 1));
band4 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band4(6:8,6:8,index), 1)), 1));
band5 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band5(6:8,6:8,index), 1)), 1));
band6 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band6(6:8,6:8,index), 1)), 1));
band7 = transpose(nanmean(squeeze(nanmean(reflectanceLandBand.Band7(6:8,6:8,index), 1)), 1));

band1Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band1(6:8,6:8,index), 1)), 1));
band2Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band2(6:8,6:8,index), 1)), 1));
band3Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band3(6:8,6:8,index), 1)), 1));
band4Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band4(6:8,6:8,index), 1)), 1));
band5Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band5(6:8,6:8,index), 1)), 1));
band6Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band6(6:8,6:8,index), 1)), 1));
band7Std = transpose(nanstd(squeeze(nanmean(reflectanceLandBand.Band7(6:8,6:8,index), 1)), 1));
%% Extract for the days that matches with PROSAIL simulation, which was based on the availability of 
% of OCO-2 SIF

% arrange to match the oco-2 time period from September 2014 to August 2015
eightDayVec = 4:8:4+45*8;   % note the assymetry - 3 days before and 4 days after. so this is not the exect mid-point of 8-day period
%eightDayVec = [eightDayVec, 365];

eodMonth = cumsum(eomday(2014,1:12));
indexDay = find(eightDayVec>eodMonth(8));  % Start of september
eightDayVec = [eightDayVec(indexDay), eightDayVec(1:(indexDay(1)-1))];  

% eight dayvec goes from day 244 (241-248) in 2014 to 236 (233-240)

% smooth to minimze noise and fill two NANs

band1 = smooth(band1(31:76,1), 3);
band2 = smooth(band2(31:76,1),3);
band3 = smooth(band3(31:76,1), 3);
band4 = smooth(band4(31:76,1), 3);
band5 = smooth(band5(31:76,1),3);
band6 = smooth(band6(31:76,1),3);
band7 = smooth(band7(31:76,1),3);

band1Std = smooth(band1Std(31:76,1), 3);
band2Std = smooth(band2Std(31:76,1),3);
band3Std = smooth(band3Std(31:76,1), 3);
band4Std = smooth(band4Std (31:76,1), 3);
band5Std = smooth(band5Std(31:76,1),3);
band6Std = smooth(band6Std(31:76,1),3);
band7Std = smooth(band7Std (31:76,1),3);

%% Collect all simulated band data

simResultDir = ['/Users/mverma/Documents/jpl/prosail/simulation-results/', N_value_dir]; 

cd (simResultDir);
fileList = dir('day_*');

load(fileList(1).name);  % to get the number of records
recordsNu = length(allResults);

simulatedBand1 = nan(recordsNu,46);
simulatedBand2 = nan(recordsNu,46);
simulatedBand3 = nan(recordsNu,46);
simulatedBand4 = nan(recordsNu,46);
simulatedBand5 = nan(recordsNu,46);
simulatedBand6 = nan(recordsNu,46);
simulatedBand7 = nan(recordsNu,46);

for n = 1:46
    
    load(fileList(n).name);
    allResults = [allResults{1:end,1}]; 
    simulatedBand1(:,n) = transpose(allResults(1,:));
    simulatedBand2(:,n) = transpose(allResults(2,:));
    simulatedBand3(:,n) = transpose(allResults(3,:));
    simulatedBand4(:,n) = transpose(allResults(4,:));
    simulatedBand5(:,n) = transpose(allResults(5,:));
    simulatedBand6(:,n) = transpose(allResults(6,:));
    simulatedBand7(:,n) = transpose(allResults(7,:));    
end;

%% Create simulated ndvi and evi - note that this is in part repetition of above, but i want to keep band and ndvi separate  

cd (simResultDir);

fileList = dir('day_*');
simulatedNdvi = nan(recordsNu, 46);
simulatedEvi = nan(recordsNu, 46);
allParameters = nan(recordsNu, 14,46);

for n = 1:46
    
    load(fileList(n).name);
    
    allResults = [allResults{1:end,1}]; 
    
    redBand = transpose(allResults(1,:));
    nirBand = transpose(allResults(2,:));
    blueBand = transpose(allResults(3,:));
    
    ndviVal = (nirBand-redBand)./(nirBand + redBand);
    eviVal = (2.5)* (nirBand-redBand)./ (nirBand + 6*redBand - 7*blueBand + 1);
    
    ndviVal = smooth(ndviVal, 3);
    eviVal = smooth(eviVal, 3);
    
    simulatedNdvi(:,n) = ndviVal;
    simulatedEvi(:,n) = eviVal;
    allParameters(:,:,n) = allCombinations;
end;

clear redBand nirBand ndviVal allCombinations allResults n fileList soilRef;

%% Calculate absolute deviation of simulated band from modis band

devBand1 = abs(simulatedBand1 - (repmat(transpose(band1(1:46)), size(simulatedBand1,1), 1)));
devBand2 = abs(simulatedBand2 - (repmat(transpose(band2(1:46)), size(simulatedBand1,1), 1)));
devBand3 = abs(simulatedBand3 - (repmat(transpose(band3(1:46)), size(simulatedBand1,1), 1)));
devBand4 = abs(simulatedBand4 - (repmat(transpose(band4(1:46)), size(simulatedBand1,1), 1)));
devBand5 = abs(simulatedBand5 - (repmat(transpose(band5(1:46)), size(simulatedBand1,1), 1)));
devBand6 = abs(simulatedBand6 - (repmat(transpose(band6(1:46)), size(simulatedBand1,1), 1)));
devBand7 = abs(simulatedBand7 - (repmat(transpose(band7(1:46)), size(simulatedBand1,1), 1)));

devNdvi = abs(simulatedNdvi - (repmat((meanNdvi(1:46)), size(simulatedBand1,1), 1)));
devEvi = abs(simulatedEvi - (repmat((meanEvi(1:46)), size(simulatedBand1,1), 1)));

% note that we do not have modis data for two 8-day periods during growing
% season

totalDev = devBand1 + devBand2 + devBand3; + devBand4  + devBand5 + devBand6 + devBand7;

b123Dev = devBand1 + devBand2 + devBand3;

[~, rowIndex1] = sort(devBand1);
[~, rowIndex2] = sort(devBand2);
[~, rowIndex3] = sort(devBand3);
[~, rowIndex4] = sort(devBand4);
[~, rowIndex5] = sort(devBand5);
[~, rowIndex6] = sort(devBand6);
[~, rowIndex7] = sort(devBand7);


[b, rowIndex] = sort(totalDev);
[b123, rowIndex123] = sort(b123Dev);

[bNdvi, rowIndexNdvi] = sort(devNdvi);
[bEvi, rowIndexEvi] = sort(devEvi);

clear c1 c2 devFun;

%% Plot and observe

eightDayVec1 = eightDayVec;
eightDayVec1(17:end) = eightDayVec1(17:end)+365;
f = @(x) num2str(x);
xLabel = arrayfun(f, eightDayVec1, 'UniformOutput', false);

sortedBand1 = simulatedBand1(rowIndex1);
figure;
boxplot(sortedBand1(1:noOfSim,:))
hold on;
plot(1:46, nanmean(sortedBand1(1:noOfSim,:)), 'r*-')
plot(band1, 'g*-')
title('Band1')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand2 = simulatedBand2(rowIndex2);
figure;
boxplot(sortedBand2(1:noOfSim,:))
hold on;
plot(1:46, nanmean(sortedBand2(1:noOfSim,:)), 'r*-')
plot(band2, 'g*-')
title('Band2')
legend('PROSAIL-Mean', 'MODIS')
ylim([0.20 0.29])
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand3 = simulatedBand3(rowIndex3);
figure; 
boxplot(sortedBand3(1:noOfSim,:))
hold on;
plot(nanmean(sortedBand3(1:noOfSim,:)), 'r*-')
plot(band3, 'g*-')
title('Band3')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand4 = simulatedBand4(rowIndex4);
figure; 
boxplot(sortedBand4(1:noOfSim,:))
hold on;
plot(nanmean(sortedBand4(1:noOfSim,:)), 'r*-')
plot(band4, 'g*-')
title('Band4')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand5 = simulatedBand5(rowIndex5);
figure; 
boxplot(sortedBand5(1:noOfSim,:))
hold on;
plot(nanmean(sortedBand5(1:noOfSim,:)), 'r*-')
plot(band5, 'g*-')
title('Band5')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand6 = simulatedBand6(rowIndex6);
figure; 
boxplot(sortedBand6(1:noOfSim,:))
hold on;
plot(nanmean(sortedBand6(1:noOfSim,:)), 'r*-')
plot(band6, 'g*-')
title('Band6')
legend('PROSAIL-Mean', 'MODIS')
ylim([0.17 0.44])
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedBand7 = simulatedBand7(rowIndex7);
figure; 
boxplot(sortedBand7(1:noOfSim,:))
hold on;
plot(nanmean(sortedBand7(1:noOfSim,:)), 'r*-')
plot(band7, 'g*-')
title('Band7')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))
    

sortedNdvi = simulatedNdvi(rowIndexNdvi);
figure; 
boxplot(sortedNdvi(1:noOfSim,:))
hold on;
plot(nanmean(sortedNdvi(1:noOfSim,:)), 'r*-')
plot(meanNdvi, 'g*-')
title('NDVI')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))

sortedEvi = simulatedEvi(rowIndexEvi);
figure;
boxplot(sortedEvi(1:noOfSim,:))
hold on;
plot(nanmean(sortedEvi(1:noOfSim,:)), 'r*-')
plot(meanEvi, 'g*-')
title('EVI')
legend('PROSAIL-Mean', 'MODIS')
set(gca, 'xtick', 1:5:46);
set(gca, 'xticklabel', xLabel(1:5:46))
set(gca, 'fontsize', 12, 'fontweight', 'bold');
xlabel('Day since Jan 1, 2014')
ylabel('EVI')
%% Parameters are in the following sequence

% N - leaf structure parameter
% Cab - Chlorophyll content
% Car - Carotenoid content
% Cbrown - Brown matter
% Cw - equivalent water thickness
% Cm - LMA
% LIDFa,LIDFb,TypeLidf - for leaf angle distribution, prescribed
% LAI - leaf area index
% hspot,tts(jj),tto,psi - are prescribed;

%% Reshape allParameters ( i do not want to use reshape)

allParametersNew = nan(recordsNu, 46, 14);

for kk=1:14
    
    temp = squeeze(allParameters(:,kk,:));
    allParametersNew (:,:,kk) = temp;
end

%% find all values that are within the margin of error of modis reflectance
% One use the total reflectance or with only ndvi or only evi or some other comination 
% Find all possible combinations that are within the margin of error of
% total including bands, ndvi, and evi (eventually I am using evi, as it
% gives best results)

% but some results are not within the margin of errors to instead I am just
% taking the 100 closest results. 


meanOptimizedParameters = nan(46, 14);
stdOptimizedParameters = nan(46, 14);

for k=1:46
    
%     tmp = find(bEvi(:,k)> 3*stdEvi(:,k), 1);
%     tmp = rowIndexNdvi(1:tmp,k);
    
    tmp = rowIndexEvi(1:noOfSim,k);
    
    temp = squeeze(allParametersNew (:,k,:));
    temp = temp (tmp,:); 
    
    meanOptimizedParameters(k,:) = nanmean(temp);
    stdOptimizedParameters (k,:) = nanstd(temp)./sqrt(noOfSim);
end    
    
    
% nRows = round(0.1*recordsNu);
% optimizedParameters = nan(nRows, 46, 14);
% 
% % Let us choose the first 3000 
% 
% for i=1:14
%     
%     temp = allParametersNew (:,:,i);
%     temp = temp (rowIndexEvi(1:nRows,:));
%     optimizedParameters(:,:,i) = temp;
% end
% 

%% Get all optimized parameters for the subseqnent analyses in SCOPE

%% Do some plotting and observe
laiVal = meanOptimizedParameters(:,10);
laiStd = stdOptimizedParameters(:,10);

chlVal = meanOptimizedParameters(:,2);
chlStd = stdOptimizedParameters(:,2);

dayVec = [eightDayVec(1:16), eightDayVec(17:end)+365];

figure;
[h, a1, a2] = plotyy(dayVec, laiVal, dayVec, chlVal);

a1.Marker = '*';
a1.MarkerSize = 12;
a1.LineWidth = 2;
h(1).YLabel.String = 'LAI';
h(1).XLabel.String = 'DAY OF YEAR';
set(h(1), 'YLim', [0 1.6]);
set(h(1), 'XLim', [dayVec(1)-10 dayVec(end)+10]);
set(h(1), 'YTick', 0:2/5:2);


a2.Marker = 'O';
a2.MarkerSize = 10;
a2.LineWidth = 1;
h(2).YLabel.String = 'CHL';
set(h(2), 'YLim', [0 55]);
set(h(2), 'XLim', [dayVec(1)-10 dayVec(end)+10]);
set(h(1), 'YTick', 0:2/5:2);
set(gca, 'xtick', 1:8:1+45*8);

    

