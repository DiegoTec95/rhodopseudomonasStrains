load('modelCGA009.mat')


model = modelCGA009;

[a,b] = exchangeSingleModel(model);

%active exchange reactions as substrates contained in variable b
%general model CGA009
%photons are controlled with 'EX_photon590_e','EX_photon630_e','EX_photon650_e','EX_photon690_e'
%CO2 uptake and excretion is regulated with the reaction 'EX_co2_e'
%mineral uptake rates can be adjusted modifying the upper or lower bound.
%However, their constraints usually remain from -1000 to 1000, unless you
%want to measure the impact of low mineral concentrations

%minerals available in the model for Rhodopseudomonas 'EX_co2_e','EX_cobalt2_e','EX_zn2_e','EX_so4_e'
%'EX_ca2_e','EX_mn2_e','EX_mg2_e','EX_cu2_e','EX_k_e','EX_fe3_e','EX_mobd_e','EX_na1_e','EX_cl_e','EX_bo3_e'


%nitrogen source
%currently the model is under anaerobic diazotrophic conditions
%nitrogen fixation n2_e
%if ammonium is required to be the unique nitrogen source -> switch EX_n2_e
%to EX_nh4_e

%carbon source for CGA009 -> EX_co2_e and EX_ac_e
%this model (CGA009) and TIE1 can perform photosynthesis



%fasta file with protein sequences of CGA009 = GCF_prot.faa


%fasta file with protein sequences of TIE1 = GCF_TIE1.faa 

%CGA669
%Rubisco knock out (CGA669) 
%stops photosynthesis 


%to switch to aerobic conditions,EX_o2_e lower bound has to be set to negative values 


%to determine the experimental uptake rates, please consider:

%for example -> acetate 

% uptake rate mmol/gDWh -> mmol/grams of Dry Weight per hour

% we need the initial and final concentration of the nutrient, in this
% case, acetate (concentration must be in mmol).

% we also need the total production of biomass = final and initial
% concentration in grams of dry weight

%finally, we need the time or the duration of the experiment, for example
%72 h

% the final operation would be (initAcetate - finalAcetate)/((final
% DW-initial DW)(time in hours))

%you can check the available nutrients that the model can use or produce
%based on the exchange reactions (look for the reactions that starts with 'EX_')


%relevant files: 

%CGA669.mat mutant model 

%modelCGA009.mat model of R. palustris CGA009

%modelTIE1.mat model of R. palustris TIE1







