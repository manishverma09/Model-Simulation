function output = prepareInputData_forScopeSimulation_atSturtPlains

% Extract input data at satellite overpass time for SCOPE simulation and
% validation

%% ######################################
% GET DATA FROM FLUX TOWER MEASUREMENTS - sw incoming, lw incoming, ea, ta,
% and GPP
% #######################################
%% Load fluxtower data
load('site_v12a_swInFparEvi.mat');
load('site_v12a.mat');

% calculate vapor pressure

eSat = 6.1094*exp(17.625*ta./(ta+243.01));
ea = eSat - vpd*10;
%% Extract at the time of satellite overpass time to run scope simulation

index1 = yearDoyTimeGpp(:,1) == 2014 | yearDoyTimeGpp(:,1) == 2015;

yearDoyTimeGpp = yearDoyTimeGpp(index1,:);

fsd = fsd(index1,:);
fld = fld(index1);
ta = ta(index1);
ea = ea(index1);

fsd(fsd<=0) = nan;
fld(fld<0) = nan;

%% Days we have SIF observations in 2014 and 2015

uniqueDays = [282,330,362,394,403,426,435,458,467,506,538,547,554,563];

%% Extract at the time of satellite overpass time

index1 = yearDoyTimeGpp(:,1) ==2015;
yearDoyTimeGpp(index1,2) = yearDoyTimeGpp(index1,2)+365;

index1 = ismember(yearDoyTimeGpp(:,2), uniqueDays);
yearDoyTimeGpp = yearDoyTimeGpp(index1,:);

fsd = fsd(index1,:);
fld = fld(index1);
ta = ta(index1);
ea = ea(index1);

%% extract data at satellite overpass time which is between 2:15 and
% 2:30 at all the stations. You can manually confirm this from sif
% data which has a field called time.

newIndex = (abs(yearDoyTimeGpp(:,3)-0.6042)<0.0001);
yearDoyTimeGpp = yearDoyTimeGpp(newIndex,:);

fsd = fsd (newIndex);
fld = fld (newIndex);
ta = ta (newIndex);
ea = ea(newIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET DATA FROM PROSAIL SIMULATION - lai, cab, and othe prospect parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prosail = load ('workspace/prosail_optimizedParameters');

cab = prosail.sifDayCab_interpolated;
lai = prosail.sifDayLai_interpolated;

% cw = interp1(prosail.dayProsail, prosail.meanOptimizedParameters(:,5), prosail.daySif);
% Note that the cw above is nearly a constant so we can take a single value
% of cw

cw = 0.010; 

cabStd = prosail.sifDayCabStd;
laiStd = prosail.sifDayLaiStd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET GEOMETRY - the only time series we have have is 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get solar zenith angle at OCO-2 overpass time for the days SIF was retrieved
daySif = prosail.daySif;
daySif(daySif>365) = daySif(daySif>365) - 365;

latitude = -17.1507; longitude = 133.3502;  % Sturt Plains
gmtDiff = 9.5;  % in hours
lstm = 15*(gmtDiff);

b = (360/365)*(daySif-81);   % in degrees
eot = 9.87*sind(2*b) - 7.53*cosd(b) - 1.5*sind(b);
tc = 4*(longitude - lstm) + eot;
lst = 14.25 + (tc/60);
hra = 15*(lst-12);


solarZenith = nan (numel(hra),1);
decAngle = nan(numel(hra),1);

for i=1:numel(hra)
    
    declination = (23.45)*(pi/180)* sin(2*pi*((284+daySif (i))/365));
    declination = declination*(180/pi);  % in degrees
    decAngle(i,1) = declination;
    
    cosZeAngle = sind(declination)*sind(latitude) + cosd(declination)*cosd(latitude)*cosd(hra(i));
    
    solarZenith(i,1) = acosd(cosZeAngle);
end

tts = solarZenith;

%% Get sensor zenith angle at OCO-2 overpass time for the days SIF was retrieved

sif = shaperead('site-specific/site.shp');

doyVal = [sif.doy];
yearVal = [sif.year];
sensorZ = [sif.sensorZeni];

index = yearVal==2015;
doyVal(index) = doyVal(index)+365;
sensorZ = grpstats(sensorZ, doyVal);

tto = sensorZ;
%% Collect all output
output.fsd = fsd;
output.fld= fld;
output.ta = ta;
output.ea = ea;
output.cab = cab;
output.lai = lai;
output.cw = cw;
output.cabStd = cabStd;
output.laiStd=laiStd;
output.tts = tts;
output.tto=tto; 