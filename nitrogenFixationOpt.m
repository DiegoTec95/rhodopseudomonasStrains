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
    sol1 = optimizeCbModel(model_temp,'max',0,0);
    max_growth = sol1.f;
    growth_rates(i) = max_growth;

    % Step 2: Constrain growth and maximize NH4 excretion
    model_temp = changeRxnBounds(model_temp, 'BIOMASS__1', growth_fraction * max_growth, 'l');
    model_temp = changeObjective(model_temp, 'EX_nh4_e');
    sol2 = optimizeCbModel(model_temp,'max',0,0);
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
