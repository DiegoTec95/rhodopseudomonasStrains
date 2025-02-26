%lets check what carbon sources are suitable for consumption using CGA009
%PM1 -> first set of carbon sources
%PM2 -> second set of carbon sources
%PM3 -> unique set of nitrogen sources
%the carbon and nitrogen sources are based on experiments performed in our
%lab. For more details, 
% please read: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011371

%initialize COBRA Toolbox and change Solver
initCobraToolbox()
changeCobraSolver('gurobi','all')

%load
dataPM1 = readcell('Supplementary Material 2.xlsx','Sheet','PM1');
dataPM2 = readcell('Supplementary Material 2.xlsx','Sheet','PM2');
dataPM3 = readcell('Supplementary Material 2.xlsx','Sheet','PM3');

reducedPM1 = dataPM1(2:80,2:3);
%load your model, for this scenario we will use CGA009
load('modelCGA009.mat')
model = modelCGA009;

%lets check the current constraints
[a,b] = exchangeSingleModel(model);


%under aerobic non phototrophic/non diazotrophic conditions
%PM1 initial carbon sources:
redPM1 = dataPM1(2:80,1:3);

carbonPM1 = redPM1(:,2:3);


for i = 1 : length(redPM1)
    %aerobic conditions => no CO2 assimilation, no light, no n2 uptake, o2
    %active, and nh4 uptake
    modelPM1 = changeRxnBounds(model,'EX_co2_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_n2_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon590_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon630_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon650_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon690_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_nh4_e',-1.5,'l'); %adjust this value to the exp measurement
    modelPM1 = changeRxnBounds(modelPM1,'EX_o2_e',-5,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',redPM1{i,3});
    checkExch = findRxnIDs(modelPM1,exchRxn);
    if checkExch ~= 0
        modelPM1 = changeRxnBounds(modelPM1,exchRxn,-3,'l'); %this value has to be adjusted to exp measurements
        growthPM1 = optimizeCbModel(modelPM1,'max',0,0); %maximize biomass production and does not allow
        %loopless solutions (check COBRA documentation)
        numPM1(i,1) = growthPM1.f;
        carbonPM1{i,3} = growthPM1.f;%positive values (above zero) represents growth of the this strain in that
        %specific condition
    else
        numPM1(i,1) = 0;
        carbonPM1{i,3} = 0;
    end
    
end

filterPM1 = find(numPM1>0);
[maxPM1,idxPM1] = max(numPM1);

compoundsPM1 = carbonPM1(filterPM1,:);
maxGrowthPM1 = carbonPM1(idxPM1,:); %considering per mole growth, L-proline offers the highest growth rate. However,
%an analysis based on g/L should be used or an approach using concentration
%per carbon molecule to have more realistic results.



%PM2 (second set of carbon sources)
redPM2 = dataPM2(2:65,1:3);

carbonPM2 = redPM2(:,2:3);


for i = 1 : length(redPM2)
    %aerobic conditions => no CO2 assimilation, no light, no n2 uptake, o2
    %active, and nh4 uptake
    modelPM2 = changeRxnBounds(model,'EX_co2_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_n2_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon590_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon630_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon650_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon690_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_nh4_e',-1.5,'l'); %adjust this value to the exp measurement
    modelPM2 = changeRxnBounds(modelPM2,'EX_o2_e',-5,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',redPM2{i,3});
    checkExch = findRxnIDs(modelPM2,exchRxn);
    if checkExch ~= 0
        modelPM2 = changeRxnBounds(modelPM2,exchRxn,-3,'l'); %this value has to be adjusted to exp measurements
        growthPM2 = optimizeCbModel(modelPM2,'max',0,0); %maximize biomass production and does not allow
        %loopless solutions (check COBRA documentation)
        numPM2(i,1) = growthPM2.f;
        carbonPM2{i,3} = growthPM2.f;%positive values (above zero) represents growth of the this strain in that
        %specific condition
    else
        numPM2(i,1) = 0;
        carbonPM2{i,3} = 0;
    end
    
end

filterPM2 = find(numPM2>0);
[maxPM2,idxPM2] = max(numPM2);

compoundsPM2 = carbonPM2(filterPM2,:);
maxGrowthPM2 = carbonPM2(idxPM2,:); %considering per mole growth, N-Acetyl-L-Glutamic Acid offers the highest growth rate. However,
%an analysis based on g/L should be used or an approach using concentration
%per carbon molecule to have more realistic results.


%lets merge both PM1 and PM2 lists

carbonSources = vertcat(compoundsPM1,compoundsPM2);



%PM3 nitrogen sources under aerobic conditions:
redPM3 = dataPM3(2:90,1:2);

nitrogenPM3 = redPM3;


for i = 1 : length(redPM3)
    %aerobic conditions => no CO2 assimilation, no light, no n2 uptake, o2
    %active, and nh4 uptake
    modelPM3 = changeRxnBounds(model,'EX_co2_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_n2_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_photon590_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_photon630_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_photon650_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_photon690_e',0,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_nh4_e',0,'l'); 
    modelPM3 = changeRxnBounds(modelPM3,'EX_o2_e',-5,'l');
    modelPM3 = changeRxnBounds(modelPM3,'EX_ac_e',-3,'l');  %adjust this value to the exp measurement
    exchRxn = strcat('EX_',redPM3{i,2});
    checkExch = findRxnIDs(modelPM3,exchRxn);
    if checkExch ~= 0
        modelPM3 = changeRxnBounds(modelPM3,exchRxn,-3,'l'); %this value has to be adjusted to exp measurements
        growthPM3 = optimizeCbModel(modelPM3,'max',0,0); %maximize biomass production and does not allow
        %loopless solutions (check COBRA documentation)
        numPM3(i,1) = growthPM3.f;
       nitrogenPM3{i,3} = growthPM3.f;%positive values (above zero) represents growth of the this strain in that
        %specific condition
    else
        numPM3(i,1) = 0;
        nitrogenPM3{i,3} = 0;
    end
    
end

filterPM3 = find(numPM3>0);
[maxPM3,idxPM3] = max(numPM3);

compoundsPM3 = nitrogenPM3(filterPM3,:);
maxGrowthPM3 = nitrogenPM3(idxPM3,:); %considering per mole growth, L-proline offers the highest growth rate. However,
%an analysis based on g/L should be used or an approach using concentration
%per carbon molecule to have more realistic results.


%DIAZOTROPHIC ANAEROBIC CONDITIONS (DARK CONDITIONS)

%PM1
carbonPM1_Diazo = redPM1(:,2:3);


for i = 1 : length(redPM1)
    %aerobic conditions => no CO2 assimilation, no light, no n2 uptake, o2
    %active, and nh4 uptake
    modelPM1 = changeRxnBounds(model,'EX_co2_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_n2_e',-0.5,'l'); %adjust this value to the exp measurement
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon590_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon630_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon650_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_photon690_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_nh4_e',0,'l'); 
    modelPM1 = changeRxnBounds(modelPM1,'EX_o2_e',0,'l');
    modelPM1 = changeRxnBounds(modelPM1,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',redPM1{i,3});
    checkExch = findRxnIDs(modelPM1,exchRxn);
    if checkExch ~= 0
        modelPM1 = changeRxnBounds(modelPM1,exchRxn,-1.5,'l'); %this value has to be adjusted to exp measurements
        growthPM1 = optimizeCbModel(modelPM1,'max',0,0); %maximize biomass production and does not allow
        %loopless solutions (check COBRA documentation)
        numPM1_Diazo(i,1) = growthPM1.f;
        carbonPM1_Diazo{i,3} = growthPM1.f;%positive values (above zero) represents growth of the this strain in that
        %specific condition
    else
        numPM1_Diazo(i,1) = 0;
        carbonPM1_Diazo{i,3} = 0;
    end
    
end

filterPM1_Diazo = find(numPM1_Diazo>0);
[maxPM1_Diazo,idxPM1_Diazo] = max(numPM1_Diazo);

compoundsPM1_Diazo = carbonPM1_Diazo(filterPM1_Diazo,:);
maxGrowthPM1_Diazo = carbonPM1_Diazo(idxPM1_Diazo,:); %considering per mole growth, D-Cellobiose offers the highest growth rate. However,
%an analysis based on g/L should be used or an approach using concentration
%per carbon molecule to have more realistic results.


%PM2 (second set of carbon sources)

carbonPM2_Diazo = redPM2(:,2:3);


for i = 1 : length(redPM2)
    %aerobic conditions => no CO2 assimilation, no light, no n2 uptake, o2
    %active, and nh4 uptake
    modelPM2 = changeRxnBounds(model,'EX_co2_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_n2_e',-0.5,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon590_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon630_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon650_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_photon690_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_nh4_e',0,'l'); %adjust this value to the exp measurement
    modelPM2 = changeRxnBounds(modelPM2,'EX_o2_e',0,'l');
    modelPM2 = changeRxnBounds(modelPM2,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',redPM2{i,3});
    checkExch = findRxnIDs(modelPM2,exchRxn);
    if checkExch ~= 0
        modelPM2 = changeRxnBounds(modelPM2,exchRxn,-1.5,'l'); %this value has to be adjusted to exp measurements
        growthPM2 = optimizeCbModel(modelPM2,'max',0,0); %maximize biomass production and does not allow
        %loopless solutions (check COBRA documentation)
        numPM2_Diazo(i,1) = growthPM2.f;
        carbonPM2_Diazo{i,3} = growthPM2.f;%positive values (above zero) represents growth of the this strain in that
        %specific condition
    else
        numPM2_Diazo(i,1) = 0;
        carbonPM2_Diazo{i,3} = 0;
    end
    
end

filterPM2_Diazo = find(numPM2_Diazo>0);
[maxPM2_Diazo,idxPM2_Diazo] = max(numPM2_Diazo);

compoundsPM2_Diazo = carbonPM2_Diazo(filterPM2_Diazo,:);
maxGrowthPM2_Diazo = carbonPM2_Diazo(idxPM2_Diazo,:);
%lets merge both PM1 and PM2 lists

carbonSources_Diazo = vertcat(compoundsPM1_Diazo,compoundsPM2_Diazo);


%Example of optimization pipelines for ammonium excretion under
%diazotrophic conditions

 % Define the weights (adjust as needed)
alpha = 0.75; % Weight for BOF -> this is the percentage of the optimum growth rate compared to the wild type
beta = 0.25;  % Weight for EX_nh4_e -> this refers to the percentage of ammonium excreted based on the total fixed nitrogen


optimized_scores = zeros(length(carbonSources_Diazo),1); % Store weighted objective values
growth_rates = zeros(length(carbonSources_Diazo),1);
nh4_excretion_rates = zeros(length(carbonSources_Diazo),1);

for i = 1:length(carbonSources_Diazo)
    % Reset carbon uptake
    model_temp = model; % Copy model for each iteration
    model_temp = changeRxnBounds(model_temp,'EX_co2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_n2_e',-0.5,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon590_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon630_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon650_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon690_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_nh4_e',0,'l'); %adjust this value to the exp measurement
    model_temp = changeRxnBounds(model_temp,'EX_o2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',carbonSources_Diazo{i,2});
    model_temp = changeRxnBounds(model_temp,exchRxn,-1.5,'l'); %adjust this value to the exp measurement
    % Set combined objective function
    model_temp = changeObjective(model_temp, {'BIOMASS__1', 'EX_nh4_e'}, [alpha, beta]);
    
    % Optimize model
    solution = optimizeCbModel(model_temp,'max',0,0);
    
    % Store results
    optimized_scores(i) = solution.f; % Combined objective function value
    
    % Get individual fluxes
    model_temp = changeObjective(model_temp, 'BIOMASS__1');
    growth_rates(i) = optimizeCbModel(model_temp,'max',0,0).f;
    
    model_temp = changeObjective(model_temp, 'EX_nh4_e');
    nh4_excretion_rates(i) = optimizeCbModel(model_temp,'max',0,0).f;
end

% Sort results based on optimized scores
[sorted_scores, sort_idx] = sort(optimized_scores, 'descend');
sorted_carbon_sources = carbonSources_Diazo(sort_idx);
sorted_growth_rates = growth_rates(sort_idx);
sorted_nh4_excretion = nh4_excretion_rates(sort_idx);

% Display results
disp('Carbon source ranking based on combined growth and NH4 excretion:');
for j = 1:length(sorted_carbon_sources)
    fprintf('%s: Score = %.3f, Growth = %.3f, NH4 Excretion = %.3f\n', ...
        sorted_carbon_sources{j}, sorted_scores(j), sorted_growth_rates(j), sorted_nh4_excretion(j));
end



%Second approach 


growth_fraction = 0.85; % Set growth to 95% of max for each carbon source
nh4_excretion_rates = zeros(length(carbonSources_Diazo),1); % Store NH4 excretion rates
growth_rates = zeros(length(carbonSources_Diazo),1);

for i = 1:length(carbonSources_Diazo)
    % Reset carbon uptake
    model_temp = model; % Copy model for each iteration
    model_temp = changeRxnBounds(model_temp,'EX_co2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_n2_e',-0.5,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon590_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon630_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon650_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon690_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_nh4_e',0,'l'); %adjust this value to the exp measurement
    model_temp = changeRxnBounds(model_temp,'EX_o2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',carbonSources_Diazo{i,2});
    model_temp = changeRxnBounds(model_temp,exchRxn,-1.5,'l'); %adjust this value to the exp measurement

    % Step 1: Maximize Growth
    model_temp = changeObjective(model_temp, 'BIOMASS__1');
    sol1 = optimizeCbModel(model_temp);
    max_growth = sol1.f;
    growth_rates(i) = max_growth;

    % Step 2: Constrain growth and maximize NH4 excretion
    model_temp = changeRxnBounds(model_temp, 'BIOMASS__1', growth_fraction * max_growth, 'l');
    model_temp = changeObjective(model_temp, 'EX_nh4_e');
    sol2 = optimizeCbModel(model_temp);
    nh4_excretion_rates(i) = sol2.f;
end

% Sort based on NH4 excretion under 95% max growth
[sorted_nh4, sort_idx] = sort(nh4_excretion_rates, 'descend');
sorted_carbon_sources = carbonSources_Diazo(sort_idx);
sorted_growth_rates = growth_rates(sort_idx);

% Display results
disp('Carbon source ranking based on NH4 excretion at 95% max growth:');
for j = 1:length(sorted_carbon_sources)
    fprintf('%s: Growth = %.3f, NH4 Excretion = %.3f\n', ...
        sorted_carbon_sources{j}, sorted_growth_rates(j), sorted_nh4(j));
end


%Third approach using FVA 

% Initialize result storage
growth_rates = zeros(length(carbonSources_Diazo),1);
min_nh4_excretion = zeros(length(carbonSources_Diazo),1);
max_nh4_excretion = zeros(length(carbonSources_Diazo),1);

growth_fraction = 0.90; % Constrain growth to 95% of max

for i = 1:length(carbonSources_Diazo)
    % Reset carbon uptake
    model_temp = model; % Copy model for each iteration
    model_temp = changeRxnBounds(model_temp,'EX_co2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_n2_e',-0.5,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon590_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon630_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon650_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_photon690_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_nh4_e',0,'l'); %adjust this value to the exp measurement
    model_temp = changeRxnBounds(model_temp,'EX_o2_e',0,'l');
    model_temp = changeRxnBounds(model_temp,'EX_ac_e',0,'l');
    exchRxn = strcat('EX_',carbonSources_Diazo{i,2});
    model_temp = changeRxnBounds(model_temp,exchRxn,-1.5,'l'); %adjust this value to the exp measurement
    % Step 1: Find maximum growth rate
    model_temp = changeObjective(model_temp, 'BIOMASS__1');
    sol1 = optimizeCbModel(model_temp,'max',0,0);
    max_growth = sol1.f;
    growth_rates(i) = max_growth;

    % Step 2: Constrain growth to 95% of max
    model_temp = changeRxnBounds(model_temp, 'BIOMASS__1', growth_fraction * max_growth, 'l');

    % Step 3: Run FVA for ammonium excretion (`EX_nh4_e`)
    [minFlux, maxFlux] = fluxVariability(model_temp, 100, 'max', 'EX_nh4_e', 0, 0 );
    
    % Store FVA results
    min_nh4_excretion(i) = minFlux;
    max_nh4_excretion(i) = maxFlux;
end

% Sort based on maximum NH4 excretion under 95% max growth
[sorted_max_nh4, sort_idx] = sort(max_nh4_excretion, 'descend');
sorted_carbon_sources = carbonSources_Diazo(sort_idx);
sorted_growth_rates = growth_rates(sort_idx);
sorted_min_nh4 = min_nh4_excretion(sort_idx);

% Display results
disp('Carbon source ranking based on NH4 excretion at 95% max growth:');
for j = 1:length(sorted_carbon_sources)
    fprintf('%s: Growth = %.3f, NH4 Excretion Range = [%.3f, %.3f]\n', ...
        sorted_carbon_sources{j}, sorted_growth_rates(j), sorted_min_nh4(j), sorted_max_nh4(j));
end


%Methods 2 and 3 usually are the most efficient ones. Method 1 lacks of
%precision.



%Transcriptomics integration with duplicates

path = '/mnt/disk8TB/CitrusClementina/rootstockData/Trimmed_Correct/counts/'; %this is the path to your transcriptomics data!
fullFileName = fullfile(path, 'normalized_counts.csv'); %normalized counts (could be TPM, CPM or DESeq)
counts = readtable(fullFileName);
counts = renamevars(counts,"Var1","Genes");
metaPath = fullfile(path,'metadataRoot.csv'); %metadata file
metadata = readtable(metaPath);
genesRNA = counts.Genes;
redCounts = counts(:,2:end);
cellCounts = table2array(redCounts);
model1 = model;
%FallUSSR -> this is the condition that you want to study or analyze
fallUSSR = [1 14]; %these are the indexes with the normalized counts that correspond to the sample you want to study
%in this specific condition, the column 1 and 14 are the columns with the
%normalized counts (if you have triplicates I can provide an updated code to handle triplicates or more replicates)
modelSubs =model1;
%Here we detect all the pathways involved in the model to cluster the
%fluxes later based on the subystems
allSubs2 = unique(modelSubs.subSystems);
distSubs = allSubs2;
idxSubs = cellfun(@(x) strmatch(x,modelSubs.subSystems,"exact"),distSubs,'UniformOutput',false);
distSubs(:,2) = cellfun(@(x) modelSubs.rxns(x),idxSubs,'UniformOutput',false);
distSubs(:,3) = cellfun(@(x) findGenesFromRxns(modelSubs,x),distSubs(:,2),'UniformOutput',false);
distSubs(:,4) = cellfun(@(x) unique(vertcat(x{:})),distSubs(:,3),'UniformOutput',false);
nonEmptyIdx = find(~cellfun(@isempty,distSubs(:,4)));
redDistSubs = distSubs(nonEmptyIdx,:);



modelfallUSSR = model1;
optfallUSSR = optimizeCbModel(modelfallUSSR,'max',0,0);
actRxnsfallUSSR = modelfallUSSR.rxns(optfallUSSR.x~=0);
actRxnsfallUSSR(:,2) = num2cell(abs(optfallUSSR.x(optfallUSSR.x~=0)));
genesActfallUSSR = findGenesFromRxns(modelfallUSSR,actRxnsfallUSSR(:,1));
uniqGenesfallUSSR = unique(vertcat(genesActfallUSSR{:}));
nonemptyfallUSSR = find(~cellfun(@isempty,genesActfallUSSR));
redRxnsfallUSSR = actRxnsfallUSSR(nonemptyfallUSSR,:);
redGenesfallUSSR = genesActfallUSSR(nonemptyfallUSSR,:);
for i = 1 : length(uniqGenesfallUSSR)
    findG = cellfun(@(x) strcmp(x,uniqGenesfallUSSR{i,1}),redGenesfallUSSR,'UniformOutput',false);
    sumG = find(cellfun(@(x) sum(x),findG)>0);
    combFlux = sum(cell2mat(actRxnsfallUSSR(sumG,2)));
    uniqGenesfallUSSR{i,2} = combFlux; 
end

modelfallUSSR = generateRules(modelfallUSSR);
[grRatiofallUSSR, grRateKOfallUSSR, grRateWTfallUSSR, hasEffectfallUSSR, delRxnsfallUSSR, fluxSolutionfallUSSR] = singleGeneDeletion(modelfallUSSR, 'FBA');
grRatiofallUSSR = fillmissing(grRatiofallUSSR, 'constant', 0);
noEffectIdxfallUSSR = find(grRatiofallUSSR>=0.95); 
downRIdxfallUSSR = find(grRatiofallUSSR<0.95 & grRatiofallUSSR>=0.01);
lethalGIdxfallUSSR = find(grRatiofallUSSR<0.01);
lethalGfallUSSR = modelfallUSSR .genes(lethalGIdxfallUSSR);
%mLethalfallUSSR = intersect(missingGenes,lethalGfallUSSR);
fallUSSRIdx = fallUSSR;
fallUSSR = cellCounts(:,fallUSSRIdx);
fallUSSR(:,3) = sum(fallUSSR,2);
fallUSSR(:,4) = mean(fallUSSR(:,1:2),2);
h = histogram(fallUSSR(:,4),10000);
h.BinLimits = [0 5000];
h.BinWidth = 5;
xlim([0 5000])
commonfallUSSR = genesRNA(idxRNA);
commonfallUSSR(:,2:3) = num2cell(fallUSSR(idxRNA,2:3));
redfallUSSR = commonfallUSSR(:,1);
for j = 2 : size(commonfallUSSR,2)
    for i = 1 : length(redfallUSSR)
        findG = cellfun(@(x) strcmp(x,commonfallUSSR{i,1}),redDistSubs(:,4),'UniformOutput',false);
        sumG = find(cellfun(@(x) sum(x),findG)>0);
        %combFlux = sum(cell2mat(actRxnsC24(sumG,2)));
        redfallUSSR{i,j} = redDistSubs(sumG,1);
    end
end
nonEmptyIdxfallUSSR = find(~cellfun(@isempty,redfallUSSR(:,3)));
redfallUSSR = redfallUSSR(nonEmptyIdxfallUSSR,:);
redDistSubsfallUSSR = redDistSubs;


totalActfallUSSR = unique(vertcat(commonfallUSSR(:,1),uniqGenesfallUSSR(:,1)));
totalActfallUSSR(:,2:4) = num2cell(zeros(length(totalActfallUSSR(:,1)),3));
[~,idx2,idx2_2] = intersect(commonfallUSSR(:,1),totalActfallUSSR(:,1));
[~,idx3,idx3_2] = intersect(commonfallUSSR(:,1),totalActfallUSSR(:,1));
%[~,idx4,idx4_2] = intersect(commonfallUSSR(:,1),totalActfallUSSR(:,1));
[~,idx5,idx5_2] = intersect(uniqGenesfallUSSR(:,1),totalActfallUSSR(:,1));
totalActfallUSSR(idx2_2,2) = commonfallUSSR(idx2,2);
totalActfallUSSR(idx3_2,3) = commonfallUSSR(idx3,3);
%totalActfallUSSR(idx4_2,4) = commonfallUSSR(idx4,4);
totalActfallUSSR(idx5_2,4) = uniqGenesfallUSSR(idx5,2);


[modelfallUSSRIr, ~, ~, ~] = convertToIrreversible (modelfallUSSR);
%[c4_1,c4_2] = exchangeSingleModel(modelCornIr4);
optfallUSSRIr = optimizeCbModel(modelfallUSSRIr,'max',0);
actRxnsfallUSSRIr = modelfallUSSRIr.rxns(optfallUSSRIr.x~=0);
actRxnsfallUSSRIr(:,2) = num2cell(abs(optfallUSSRIr.x(optfallUSSRIr.x~=0)));
genesActfallUSSR2 = findGenesFromRxns(modelfallUSSRIr,actRxnsfallUSSRIr(:,1));
uniqGenesfallUSSR2 = unique(vertcat(genesActfallUSSR2{:}));
nonemptyfallUSSR2 = find(~cellfun(@isempty,genesActfallUSSR2));
redRxnsfallUSSR2 = actRxnsfallUSSRIr(nonemptyfallUSSR2,:);
redGenesfallUSSR2 = genesActfallUSSR2(nonemptyfallUSSR2,:);
for i = 1 : length(uniqGenesfallUSSR2)
    findG = cellfun(@(x) strcmp(x,uniqGenesfallUSSR2{i,1}),redGenesfallUSSR2,'UniformOutput',false);
    sumG = find(cellfun(@(x) sum(x),findG)>0);
    combFlux = sum(cell2mat(actRxnsfallUSSRIr(sumG,2)));
    uniqGenesfallUSSR2{i,2} = combFlux; 
end
totalActfallUSSR(:,5) = num2cell(zeros(length(totalActfallUSSR(:,1)),1));
[~,idx6,idx6_2] = intersect(uniqGenesfallUSSR2(:,1),totalActfallUSSR(:,1));
totalActfallUSSR(idx6_2,5) = uniqGenesfallUSSR2(idx6,2);


clear findfallUSSR

matrixfallUSSR = cell2mat(totalActfallUSSR(:,[2 3 5]));
findfallUSSR{1,1} = find(matrixfallUSSR(:,1)>0 & matrixfallUSSR(:,3)==0);
findfallUSSR{2,1} = find(matrixfallUSSR(:,2)>0 & matrixfallUSSR(:,3)==0);
%findfallUSSR{3,1} = find(matrixfallUSSR(:,3)>0 & matrixfallUSSR(:,4)==0);
clear genefallUSSR

genefallUSSR{1,1} = totalActfallUSSR(findfallUSSR{1,1},1);
genefallUSSR{2,1} = totalActfallUSSR(findfallUSSR{2,1},1);
%genefallUSSR{3,1} = totalActfallUSSR(findfallUSSR{3,1},1);

for k = 1 : length(genefallUSSR)
    k
    testedRxnsfallUSSR = {};
    modelRNAfallUSSR = modelfallUSSRIr; cont = 1;
    geneList = genefallUSSR{k,1};
    clear rxnsMod
    for i = 1 : length(geneList)
        findG = strfind(modelfallUSSRIr.grRules,geneList{i,1});
        idxCheck = find(~cellfun(@isempty,findG));
        rxnsId = modelfallUSSRIr.rxns(idxCheck);
        compRxns = setdiff(rxnsId,testedRxnsfallUSSR);
        i
        if ~isempty(compRxns)
            for j = 1 : length(compRxns)
                modelTest = changeRxnBounds(modelRNAfallUSSR,compRxns{j,1},0.01,'l');
                testGrowth = optimizeCbModel(modelTest,'max',0,0);
                adjGrowth = 0.9 * 0.3 * optfallUSSRIr.f;
                if testGrowth.f> adjGrowth
                    modelRNAfallUSSR = modelTest;
                    testedRxnsfallUSSR = unique(vertcat(testedRxnsfallUSSR,compRxns));
                    rxnsMod{cont,1} = compRxns{j,1};
                    cont = cont + 1;
                    break
                    
                end
            end
        
        end    
        % sumG = find(cellfun(@(x) sum(x),findG)>0); normVal = totalActC4{i,2};
        % uniqGenesC4{i,2} = combFlux;%modify this part  
    end
    modelsRNAfallUSSR{k,1} = modelRNAfallUSSR; 
    optfallUSSRIr2 = optimizeCbModel(modelRNAfallUSSR,'max',0,0);
    %[GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(modelRNAfallUSSR)
    actRxnsfallUSSRIr2 = modelRNAfallUSSR.rxns(optfallUSSRIr2.x~=0);
    actRxnsfallUSSRIr2(:,2) = num2cell(abs(optfallUSSRIr2.x(optfallUSSRIr2.x~=0)));
    actRxnsfallUSSRIr2(:,3) = modelRNAfallUSSR.subSystems(findRxnIDs(modelRNAfallUSSR,actRxnsfallUSSRIr2(:,1)));
    subsIrfallUSSR2 = unique(actRxnsfallUSSRIr2(:,3));
    for i = 1 : length(subsIrfallUSSR2)
        xData = strcmp(subsIrfallUSSR2{i,1},actRxnsfallUSSRIr2(:,3));
        sumData = sum(cell2mat(actRxnsfallUSSRIr2(xData,2)));
        subsIrfallUSSR2{i,2} = sumData;
    end
    subsIrfallUSSR2 = subsIrfallUSSR2(~contains(subsIrfallUSSR2(:,1),'Extracellular'),:);
    subsfallUSSR{k,1} = subsIrfallUSSR2;
end
actRxnsfallUSSRIr(:,3) = modelfallUSSRIr.subSystems(findRxnIDs(modelfallUSSRIr,actRxnsfallUSSRIr(:,1)));
subsIrfallUSSR = unique(actRxnsfallUSSRIr(:,3));
for i = 1 : length(subsIrfallUSSR)
    xData = strcmp(subsIrfallUSSR{i,1},actRxnsfallUSSRIr(:,3));
    sumData = sum(cell2mat(actRxnsfallUSSRIr(xData,2)));
    subsIrfallUSSR{i,2} = sumData;
end
subsIrfallUSSR = subsIrfallUSSR(~contains(subsIrfallUSSR(:,1),'Extracellular'),:);



wholeSubsfallUSSR = unique(vertcat(subsIrfallUSSR(:,1),subsfallUSSR{1,1}(:,1),subsfallUSSR{2,1}(:,1)));
wholeSubsfallUSSR(:,2:4) = num2cell(zeros(length(wholeSubsfallUSSR(:,1)),3));
[~,idxT1,idxT2] = intersect(subsfallUSSR{1,1}(:,1),wholeSubsfallUSSR(:,1));
wholeSubsfallUSSR(idxT2,2) = subsfallUSSR{1,1}(idxT1,2);
[~,idxT1,idxT2] = intersect(subsfallUSSR{2,1}(:,1),wholeSubsfallUSSR(:,1));
wholeSubsfallUSSR(idxT2,3) = subsfallUSSR{2,1}(idxT1,2);
%[~,idxT1,idxT2] = intersect(subsfallUSSR{3,1}(:,1),wholeSubsfallUSSR(:,1));
%wholeSubsfallUSSR(idxT2,4) = subsfallUSSR{3,1}(idxT1,2);

[~,idxT1,idxT2] = intersect(subsIrfallUSSR(:,1),wholeSubsfallUSSR(:,1));
wholeSubsfallUSSR(idxT2,4) = subsIrfallUSSR(idxT1,2);
wholefallUSSR = cell2mat(wholeSubsfallUSSR(:,2:end));
sumSubs = sum(wholefallUSSR,1);
wholefallUSSR_2 = wholefallUSSR./sumSubs*100;




