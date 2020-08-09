function allResults = runProsail_sturtPlain
% run prosail at sturt plain site and estimate hemispherical directional
% reflectance for sevean modis bands

%% Get leaf angle distribution - sturt plains is a grassland site
% so assume a erectophile distribution
TypeLidf=1;
% if 2-parameters LIDF: TypeLidf=1
if (TypeLidf==1)
    % LIDFa LIDF parameter a, which controls the average leaf slope
    % LIDFb LIDF parameter b, which controls the distribution's bimodality
    %	LIDF type 		a 		 b
    %	Planophile 		1		 0
    %	Erectophile    -1	 	 0
    %	Plagiophile 	0		-1
    %	Extremophile 	0		 1
    %	Spherical 	   -0.35 	-0.15
    %	Uniform 0 0
    % 	requirement: |LIDFa| + |LIDFb| < 1
    LIDFa	=	-1;
    LIDFb	=	0;
end

%% Get geometry - we want to invert for reflectance in each 8-day period

% so at the center of each 8-day period between june and july marked

latitude = -17.1507;   % Sturt Plains
hourAngle = 0;  % MCD43 gives nadir reflectance at local solar noon

eightDayVec = 4:8:4+45*8;   % note the assymetry - 3 days before and 4 days after. so this is not the exect mid-point of 8-day period
%eightDayVec = [eightDayVec, 365];

solarZenith = nan (numel(eightDayVec),1);
decAngle = nan(numel(eightDayVec),1);

for i=1:numel(eightDayVec)
    
    declination = (23.45)*(pi/180)* sin(2*pi*((284+eightDayVec(i))/365));
    declination = declination*(180/pi);  % in degrees
    decAngle(i,1) = declination;
    
    cosZeAngle = sind(declination)*sind(latitude) + cosd(declination)*cosd(latitude)*cosd(hourAngle);
    
    solarZenith(i,1) = acosd(cosZeAngle);
end

% arrange to match the oco-2 time period from September 2014 to August 2015

eodMonth = cumsum(eomday(2014,1:12));
indexDay = find(eightDayVec>eodMonth(8));  % Start of september

eightDayVec = [eightDayVec(indexDay), eightDayVec(1:(indexDay(1)-1))];
solarZenith = [solarZenith(indexDay); solarZenith(1:(indexDay(1)-1))];

%% Get soil moisture from ozflux tower data to prescribe psoil and use mean of 8-day interval

smLai = load('metAndRsData/SturtPlains_v12a_swInFparEvi');
sm = smLai.sm3cm;     
plot(sm, '*-')   % take note of maximum values
laiVal = smLai.lai;
eviVal = smLai.evi;

doy = load('ozflux/gppData/SturtPlains_v12a');
doy = doy.yearDoyTimeGpp;
doy = doy(:,1:2);

% extract to match eightDayVec above, do it for 2014 and 2015 separately and
% then combine

% For 2014
index = doy(:,1) ==2014;
doyYear1 = doy(index,:);
smYear1 = sm(index);
laiValYear1 = laiVal(index);
eviValYear1 = eviVal(index);


index1 = find(doyYear1(:,2)==eightDayVec(1)-3, 1);  % eightDayVec values are middle of eight day period, so subtract 3
smYear1 = smYear1(index1:end);
laiValYear1 = laiValYear1(index1:end);
eviValYear1 = eviValYear1(index1:end);
doyYear1 = doyYear1(index1:end,:);

% For 2015
index = doy(:,1) ==2015;
doyYear2 = doy(index,:);
smYear2 = sm(index);
laiValYear2 = laiVal(index);
eviValYear2 = eviVal(index);


index2 = find(doyYear2(:,2)==eightDayVec(end)+4, 1, 'last');  % eightDayVec values are middle of eight day period, so subtract 3
smYear2 = smYear2(1:index2);
laiValYear2 = laiValYear2(1:index2);
eviValYear2 = eviValYear2(1:index2);
doyYear2 = doyYear2(1:index2,:);

% Combine

sm = cat(1, smYear1(:), smYear2(:));
laiVal = cat(1, laiValYear1(:), laiValYear2(:));
eviVal = cat(1, eviValYear1(:), eviValYear2(:));
doyYear = cat(1, doyYear1, doyYear2);


% extract to match with 8-day data of MODIS
meanSm = nan(46,1);
meanLai = nan(46,1);
meanEvi = nan(46,1);

newEightDayVec = 1:8:1+45*8;
%newEightDayVec = [newEightDayVec, 365];

newEightDayVec = [newEightDayVec(indexDay), newEightDayVec(1:(indexDay(1)-1))];
newEightDayVec = [newEightDayVec(:);newEightDayVec(end)+8];

k=1;
for i=1:46    %[1:16, 18:46]
    
    index = doyYear(:,2)>= newEightDayVec(i) & doyYear(:,2)< newEightDayVec(i+1);
    
    meanSm(k,1) = nanmean(sm(index));
    meanLai(k,1) = nanmean(laiVal(index));
    meanEvi(k,1) = nanmean(eviVal(index));
    k=k+1;
end

%% Now create psoil - note that the soil has nearly 45% clay.

maxLimSoil = 0.3;  % Plot meanSm and see the max where we assume that soils are fully wet

psoil = nan(46,1);

psoil(meanSm<=0) = 1;
psoil(meanSm>=maxLimSoil) = 0;

psoil(meanSm>0 & meanSm<maxLimSoil) = 1 - ((meanSm(meanSm>0 & meanSm<maxLimSoil))/(maxLimSoil-0.01));

psoil(psoil>1) = 0;

%% Soil Reflectance Properties - estimate soil reflectance at each day of simulation.
% we are using soil moisture data at 3 cm from tower dataset
% to specify psoil, which measues degree of dryness.

data=dataSpec_P5B;

Rsoil1=data(:,10);   % dry soil spectra
Rsoil2=data(:,11);   % wet soil spectra

rsoil0 = nan(numel(Rsoil1),numel(psoil));

for ii=1:numel(psoil)
    rsoil0 (:,ii)=psoil(ii).*Rsoil1+(1-psoil(ii))*Rsoil2;    
end

Es=data(:,8);Ed=data(:,9);
rd=pi/180;  % for diffuse radiation calculation used later
%% PROSPECT parameters (Darvishazedh et al., 2008, RSE for a grassland site)
% Also see Feret et al., 2008, RSE, especially the min, max, and mean from four different datasets
% in Table 1.
%LEAF CHEM & STR PROPERTIES	%%

Cab		=	0:4:80;		% chlorophyll content (µg.cm-2)
Car		=	0;            % 0:3:30;			% carotenoid content (µg.cm-2)
Cbrown	=	0.0;	% brown pigment content (arbitrary units)
Cw		=	0:0.01:0.02;	% EWT (cm)
Cm		=	0:0.0025:0.01;	% LMA (g.cm-2)
N		=	1.6;	% structure coefficient run for 1.4, 1.5, and 1.6 it is a monocot leaf so this paprameter should be on the lower side
%% Specify SAIL parameters
hspot	=	0.01;
tts = solarZenith;
tto		=	0.;		% observer zenith angle (°)
psi		=	0.;

LAI = 0:0.1:2;
%% Create all possible combination of input parameters
% argument sequence: N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0

for jj = 1:46    
    sets = {N, Cab, Car, Cbrown, Cw, Cm, LIDFa,LIDFb,TypeLidf,LAI,hspot,tts(jj),tto,psi};
    [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14] = ndgrid (sets{:});
    allCombinations = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:) x10(:) x11(:) x12(:) x13(:) x14(:)];
    
    soilRef = rsoil0(:,jj);
    % %% Simulate
    %
    % bandLowerL = [620 841 459 545 1230 1628 2105];
    % bandUpperL = [670 876 479 565 1250 1652 2155];
    %
    % bandIndex = nan(7,2);
    %
    % for k=1:7
    %
    %         bandIndex(k,1) = find(data(:,1)==bandLowerL(k));
    %         bandIndex(k,2) = find(data(:,1)==bandUpperL(k));
    % end;
    
    %allResults = cell(numel(allCombinations), 2);
    
    %% Run parallely
    
    clear allResults;
    
    parfor i=1:size(allCombinations,1)
        
        % rdot: hemispherical-directional reflectance factor in viewing direction
        % rsot: bi-directional reflectance factor
        % rsdt: directional-hemispherical reflectance factor for solar incident flux
        % rddt: bi-hemispherical reflectance factor
        
        inputArg = allCombinations(i,:);
        
        [rdot,rsot,~,~]=PRO4SAIL(inputArg(1), inputArg(2), inputArg(3), inputArg(4), inputArg(5), inputArg(6),...
            inputArg(7), inputArg(8), inputArg(9), inputArg(10), inputArg(11), inputArg(12),...
            inputArg(13), inputArg(14), soilRef);        
        
        %    [rdot,rsot,~,~]=PRO4SAIL(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	direct / diffuse light	%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the direct and diffuse light are taken into account as proposed by:
        % Francois et al. (2002) Conversion of 400–1100 nm vegetation albedo
        % measurements into total shortwave broadband albedo using a canopy
        % radiative transfer model, Agronomie
        % Es = direct
        % Ed = diffuse
        
%         skyl	=	0.847- 1.61*sin((90-inputArg(12))*rd)+ 1.04*sin((90-inputArg(12))*rd)*sin((90-inputArg(12))*rd); % % diffuse radiation
%         PARdiro	=	(1-skyl)*Es;
%         PARdifo	=	(skyl)*Ed;
%         
%         % resv : directional reflectance
%         resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);        
%         
%         resvNew = [nanmean(resv(221:271)); nanmean(resv(442:477));...
%             nanmean(resv(60:80)); nanmean(resv(146:166));...
%             nanmean(resv(831:851)); nanmean(resv(1229:1253));...
%             nanmean(resv(1706:1756))];
        
        resvNew = [nanmean(rdot(221:271)); nanmean(rdot(442:477));...
            nanmean(rdot(60:80)); nanmean(rdot(146:166));...
            nanmean(rdot(831:851)); nanmean(rdot(1229:1253));...
            nanmean(rdot(1706:1756))]
        
        allResults{i,1} = resvNew;
    end
    
    saveFileName = ['prosail/simulation-results/day_', num2str(jj)];
    save (saveFileName, 'allResults', 'allCombinations', 'soilRef');
    
end
poolobj = gcp('nocreate');
delete(poolobj);
