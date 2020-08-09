%% Get data from SCOPE simulation for different Vcmax values and 
% collect them in matrix for upper and lower values 
% for Figure

cd pySCOPE/pyscope/output;


dirList = {'vcmax10', 'vcmax30', 'vcmax60', 'vcmax90', 'vcmax120', 'vcmax150', 'vcmax180'};

result = struct('noOfPoint', nan(numel(dirList),1), 'sifR', nan(numel(dirList),1),...
    'sifRmse', nan(numel(dirList),1), 'sifGppR', nan(numel(dirList),1), 'gppR', nan(numel(dirList),1),...
    'gppRmse', nan(numel(dirList),1));

%% Get upper limit
scopeSimulatedSifU = nan(14, numel(dirList));
scopeSimulatedGppU = nan(14, numel(dirList));


for i=1:numel(dirList)
    
    cd (dirList{i});
    
    %% Get fluorescence
    
    fileName = [pwd, '/upper/dat/fluorescence.dat'];
    
    fid = fopen(fileName);   % 640 to 850 nm with 1nm step
    
    fgetl(fid); fgetl(fid);
    
    flu = nan(14, 211);
    k=1;
    while ~feof(fid)
        
        l = fgetl(fid);
        l = sscanf(l, '%f');
        flu(k,:) = l';
        k=k+1;
    end
    
    fclose(fid);
    
    % Get SIF at 757 and 771 nm
    
    wvVec = 640:1:850;
    
    index1 = wvVec == 757;
    index2 = wvVec == 771;
    
    scopeSif = [flu(:,index1), flu(:,index2)];
    meanSif  = (1/2)*(scopeSif(:,1) + (scopeSif(:,2)*1.4));
    meanSif = [meanSif(12:14); meanSif(1:11)];
    
    scopeSimulatedSifU(:,i) = meanSif;
    
    %% Get gpp, apar_cab
    
    fileName = [pwd, '/upper/dat/fluxes.dat'];
    
    fid = fopen(fileName);   % this file has fluxes
    
    varName = fgetl(fid);
    varName = strsplit(varName, ' ');
    
    varUnit = fgetl(fid);
    varUnit = strsplit(varUnit, ' ');
    
    fluxOutput = nan(size(flu,1), numel(varName));
    k=1;
    while ~feof(fid)
        
        l = fgetl(fid);
        l = sscanf(l, '%f');
        fluxOutput(k,:) = l';
        k=k+1;
    end;
    
    fclose(fid);
    
    %% Plot different variables
    
    % cahnge the variables below approrpiately Actot is netphotosynthesis,
    % aPAR_Cab is par absorbed by chlorophyll
    
    y1Var = strcmp(varName, 'Actot');
    y1Vec = fluxOutput(:,y1Var);
    y1Vec = [y1Vec(12:14); y1Vec(1:11)];
    
    scopeSimulatedGppU(:,i) = y1Vec;
    cd ..;
end


%% Get lower limit
scopeSimulatedSifL = nan(14, numel(dirList));
scopeSimulatedGppL = nan(14, numel(dirList));


for i=1:numel(dirList)
    
    cd (dirList{i});
    
    %% Get fluorescence
    
    fileName = [pwd, '/lower/dat/fluorescence.dat'];
    
    fid = fopen(fileName);   % 640 to 850 nm with 1nm step  % 640 to 850 nm with 1nm step
    
    fgetl(fid); fgetl(fid);
    
    flu = nan(14, 211);
    k=1;
    while ~feof(fid)
        
        l = fgetl(fid);
        l = sscanf(l, '%f');
        flu(k,:) = l';
        k=k+1;
    end
    
    fclose(fid);
    
    % Get SIF at 757 and 771 nm
    
    wvVec = 640:1:850;
    
    index1 = wvVec == 757;
    index2 = wvVec == 771;
    
    scopeSif = [flu(:,index1), flu(:,index2)];
    meanSif  = (1/2)*(scopeSif(:,1) + (scopeSif(:,2)*1.4));
    meanSif = [meanSif(12:14); meanSif(1:11)];
    
    scopeSimulatedSifL(:,i) = meanSif;
    
    %% Get gpp, apar_cab
    
    fileName = [pwd, '/lower/dat/fluxes.dat'];
    
    fid = fopen(fileName);   % this file has fluxes
    
    varName = fgetl(fid);
    varName = strsplit(varName, ' ');
    
    varUnit = fgetl(fid);
    varUnit = strsplit(varUnit, ' ');
    
    fluxOutput = nan(size(flu,1), numel(varName));
    k=1;
    while ~feof(fid)
        
        l = fgetl(fid);
        l = sscanf(l, '%f');
        fluxOutput(k,:) = l';
        k=k+1;
    end
    
    fclose(fid);
    
    %% Plot different variables
    
    % cahnge the variables below approrpiately Actot is netphotosynthesis,
    % aPAR_Cab is par absorbed by chlorophyll
    
    y1Var = strcmp(varName, 'Actot');
    y1Vec = fluxOutput(:,y1Var);
    y1Vec = [y1Vec(12:14); y1Vec(1:11)];
    
    scopeSimulatedGppL(:,i) = y1Vec;
    cd ..;
end