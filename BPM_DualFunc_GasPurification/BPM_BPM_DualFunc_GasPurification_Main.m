%% General Code Description of BPM of Dual-Functional Gas Purification Processes and Author Introduction
%======================================================================================================================
%> @ Code Description:
%> @ File        BPM_BPM_DualFunc_GasPurification_Main.m
%> @ Brief       Advanced parameter identification module based on Bayesian inference for dual-functional gas purification processes.
%> @ Details     BPM_BPM_DualFunc_GasPurification_Main serves as the main program for Bayesian parameter estimation framework of dual-functional gas purification models,
%>               encompassing a holistic approach that integrates experimental data processing with Bayesian inference.
%>
%>               The process is characterized by the following steps:
%>               (1) Parameter Establishment: Initiation by defining parameters based on physical attributes of adsorbate and dual-functional adsorbent.
%>               (2) Bayesian Inference Optimization: Alignment of Bayesian inference settings with material properties and experimental conditions,
%>                   including optimization of steps, samples, and calibration of model predictions against actual measurements.
%>               (3) Iterative Refinement and Analysis: Ensures efficient and effective Bayesian inference iterations, leading to a comprehensive
%>                   analysis where estimated parameters from the model are juxtaposed with actual experimental data, often visualized through graphs.
%>               (4) Module for Bayesian Parameter Estimation: Facilitates the entire process of parameter estimation,
%>                   from setting initial values and options for Bayesian inference to numerically solving the model with estimated parameters
%>                   and optimizing these parameters through a rigorous Bayesian analysis.
%>
%> @ Author Introduction
%>              - Yesol Hyun (yesol2@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>
%> @ Date        September 2025
%> @ Version     1.0
%> @ Copyright   Copyright (c) 2025-2025 Yesol Hyun and Jung-Il Choi, Yonsei University, All rights reserved.
%> @ License     This project is released under the terms of the MIT License (see LICENSE file).
%======================================================================================================================
close all;clear all; clc; warning off;
rng(100,'twister')
uqlab
%% Experimental Data Preparation and Setup of Property Setup
%==========================================================================================================================================================================================================================================
%> @details     This section of BPM_BPM_DualFunc_GasPurification_Main.m encompasses the preparation of experimental data and the initial setup for Bayesian analysis in dual-functional gas purification models. The process includes:

%>              Data Loading and Processing:
%>              - Loading and calculating properties such as flow rates, bed length, particle diameter, bed volume, bulk density, and superficial velocity.
%>              - Reading and processing experimental data from CSV files, including time and concentration data.

%>              This process ensures that all the necessary steps are taken to accurately prepare the experimental data.
%===========================================================================================================================================================================================================================================
% Load measured data
cFeed_ppm= [1980 1065 570 162]; % Feed concentration[ppm]
convert= (34.08/22.4)*(273/(273+20))*(1013/1013); % Conversion
cFeed_tmp = cFeed_ppm * convert;
cFeed = cFeed_tmp * 1e-3 /34.08; % Feed concentration [mol/m^3]
NumCurve = length(cFeed_ppm); % Number of experimental data curves
L= 0.1; % Bed length [m]
Q= 1.83E-5; % Flow rates[m^3/s]
Diam= 0.03; % Diameter of the column [m]
A= (Diam/2)^2*pi;% Cross-sectional area of the column [m^2]
epsilon = 0.808; % Bed porosity

V_ads = A .* L; % Volume of the adsorbent [m^3]
rho_b = 481; % Bulk density of the adsorbent [kg/m^3]
v = Q / A; % Superficial velocity [m/s]
dp =  0.002; % Particle diameter [m]
mu = 1.32e-5; % Dynamic viscosity of the gas [Pa.s]
rho_g = 1.1012; % Density of the gas [kg/m^3]
rho = 2505; % Pellet density [kg/m^3]

dz= 0.01/10; % Grid size for the spatial discretization [m]
dt= 120; % Time step [s]
ngrid = L / dz + 1; % Number of grid points [m]
us = v; % Superficial velocity alias [m/s]

ModelNum =1;
if ModelNum ==1
    ModelName = 'Langmuir';
end
%% Step 1 of Two-Step Bayesian Framework
%% Initial Preparation for Parameter Estimation
%==========================================================================================================================================================================================================================================
%> @details     This section of BPM_BPM_DualFunc_GasPurification_Main.m encompasses the preparation of experimental data and the initial setup for Bayesian analysis in dual-functional gas purification models. The process includes:

%>              (1) Data Analysis and Preparation for Bayesian Inference:
%>                  - Performing linear interpolation, numerical integration, and normalization of data.
%>                  - Preparing data for Bayesian inference by interpolating and smoothing concentration data.
%>                  - Compiling differential concentration rate data and normalized time-to-breakthrough data.

%>              (2) Bayesian Analysis Setup:
%>                  - Defining the forward model for the Bayesian analysis using a custom function.
%>                  - Setting prior distributions for various parameters based on expected ranges.
%>                  - Configuring the Bayesian inversion solver, including the selection of a sampler and defining the number of chains and steps.
%>                  - Preparing the discrepancy model and options for the Bayesian analysis.

%>              This comprehensive approach ensures that all the necessary steps are taken to accurately prepare the experimental data and set up the Bayesian analysis for effective parameter identification in dual-functional gas purification models
%===========================================================================================================================================================================================================================================================
disp('Step 1 of 2-step BPM for Dual-Functional Gas Purification')
% Processing each experimental data curve
Ndata = 30; % Number of data points for Bayesian analysis
for i = 1:NumCurve
    % Load experimental data from CSV files
    filename = sprintf('Ortiz_error_%d.csv',cFeed_ppm(i));
    exp_inf= csvread(filename,1,0);
    
    % Extract and preprocess time and concentration data
    t_exp_tmp = exp_inf(:, 1);
    idx_tmax = find(t_exp_tmp == max(t_exp_tmp), 1, 'first');
    t_data_tmp = t_exp_tmp(1:idx_tmax) * 60; % Convert time to seconds
    c_tmp = exp_inf(1:idx_tmax, 2);
    c_tmp_convert = c_tmp * convert;
    c_tmp_convert = c_tmp_convert * 1e-3 /34.08;  %[mol/m^3]
    % Interpolate concentration data for uniform time distribution
    t_data = linspace(t_data_tmp(1), max(t_data_tmp), 20 * Ndata);
    c_data = pchip(t_data_tmp, c_tmp_convert/cFeed(i), t_data);
    
    % Further data processing for Bayesian inference
    idx1 = find(c_data==0,1,'last');
    idx1= idx1+1;
    
    dt_tmp = t_data(idx1+1:end)-t_data(idx1:end-1);
    idx3=find(c_data <= 0.4,1,'last');
    Max_C(i) = c_data(idx3);
    
    t_interp = linspace(t_data(idx1),t_data(idx3+1),Ndata);% Time vector
    c_interp = interp1(t_data(idx1:idx3+1),c_data(idx1:idx3+1)*cFeed(i),t_interp, 'linear');
    tf(i)= max(t_interp)*1.2;
    
    % Numerical differentiation for rate of change in concentration
    dt_tmp = t_interp(2) - t_interp(1);
    dc_dt_new = diff(c_interp) / dt_tmp;
    max_dcdt(i) = max(dc_dt_new);
    data_dcdt = dc_dt_new/ max_dcdt(i);
    
    % Determining specific times for a given concentration level
    specific_c = 0.01;
    specific_t = pchip(c_data(idx1:length(c_data)), t_data(idx1:length(c_data)) / 60, specific_c);
    data_tbr_tmp(i) = specific_t;
    
    % Compiling all processed data for Bayesian analysis
    if i == 1
        dcdt_tmp = data_dcdt;
        t_exp = t_interp;
        c_exp = c_interp;
    else
        dcdt_tmp = [dcdt_tmp, data_dcdt];
        t_exp = [t_exp, t_interp];
        c_exp = [c_exp, c_interp];
    end
end

% Finalizing the output data for Bayesian analysis : Normalize data
max_tbr = max(data_tbr_tmp);
data_tbr = data_tbr_tmp / max(data_tbr_tmp);
myData.y = [((Ndata - 1) / Ndata) * dcdt_tmp, (1 / Ndata) * data_tbr];
myData.Name = 'Exp data';

% Setting up the forward model for Bayesian analysis
ModelOpts_Step1.mHandle = @(par) uq_DualFunc_GasPurification_Step1(par,t_exp,c_exp,dt,ngrid,dz,cFeed,epsilon,tf,v,ModelNum,rho,rho_b,rho_g,dp,mu,NumCurve,max_dcdt,max_tbr,Ndata,Max_C);

ModelOpts_Step1.isVectorized = true;
myForwardModel_Step1 = uq_createModel(ModelOpts_Step1);
BayesOpts_Step1.ForwardModel = myForwardModel_Step1;

% Defining prior distributions for Bayesian parameters
PriorOpts_Step1.Marginals(1).Name = 'KL';
PriorOpts_Step1.Marginals(1).Type = 'Uniform';
PriorOpts_Step1.Marginals(1).Parameters = [1 4];

PriorOpts_Step1.Marginals(2).Name = 'qm';
PriorOpts_Step1.Marginals(2).Type = 'Uniform';
PriorOpts_Step1.Marginals(2).Parameters = [-1 -0.1];

PriorOpts_Step1.Marginals(3).Name = 'Dm';
PriorOpts_Step1.Marginals(3).Type = 'Uniform';
PriorOpts_Step1.Marginals(3).Parameters = [-7 -5];

PriorOpts_Step1.Marginals(4).Name = 'De';
PriorOpts_Step1.Marginals(4).Type = 'Uniform';
PriorOpts_Step1.Marginals(4).Parameters = [-7 -5];

PriorOpts_Step1.Marginals(5).Name = 'kc';
PriorOpts_Step1.Marginals(5).Type = 'Uniform';
PriorOpts_Step1.Marginals(5).Parameters = [-6 -3];

PriorOpts_Step1.Marginals(6).Name = 'kd';
PriorOpts_Step1.Marginals(6).Type = 'Uniform';
PriorOpts_Step1.Marginals(6).Parameters = [-6 -3];

% Creating the prior distribution as an input object
myPriorDist_Step1 = uq_createInput(PriorOpts_Step1);
BayesOpts_Step1.Prior = myPriorDist_Step1;

% Setting up the discrepancy model for the Bayesian analysis
SigmaOpts_Step1.Marginals.Name = 'sigma2';
SigmaOpts_Step1.Marginals.Type = 'Uniform';
SigmaOpts_Step1.Marginals.Parameters = [0 (0.2)^2];
mySigmaDist = uq_createInput(SigmaOpts_Step1);
DiscrepancyOptsUnknownDisc_Step1.Type = 'Gaussian';
DiscrepancyOptsUnknownDisc_Step1.Prior = mySigmaDist;
BayesOpts_Step1.Discrepancy = DiscrepancyOptsUnknownDisc_Step1;

% Configuring the Bayesian inversion solver and sampler
Solver_Step1.Type = 'MCMC';
Solver_Step1.MCMC.Sampler = 'AIES';
Solver_Step1.MCMC.NChains = 100;
Solver_Step1.MCMC.Visualize.Interval = 10;
Solver_Step1.MCMC.Steps =500;
Solver_Step1.MCMC.Visualize.Parameters = [1 2 3 4 5 6];

% Finalizing Bayesian inversion setup
BayesOpts_Step1.Type = 'Inversion';
BayesOpts_Step1.Data = myData;
BayesOpts_Step1.Solver = Solver_Step1;
%% Execution of Bayesian Analysis and Forward Modeling Using IUQ Results
%======================================================================================================================
%> @details     This section performs the Bayesian analysis for gas adsorption parameter identification and conducts forward modeling using the inferred parameters.
%>              The process includes:

%>              (1) Bayesian Analysis Execution:
%>                  - Perform Bayesian analysis using the specified options.
%>
%>              (2) Post-Processing of Bayesian Results:
%>                  - Filter out MCMC chains with low acceptance rates.
%>                  - Display and analyze the acceptance rates of the chains.
%>
%>              (3) Results Reporting and Visualization:
%>                  - Print a detailed report of the Bayesian analysis results.
%>                  - Display results including prior and posterior distributions.
%>
%>              (4) Forward Modeling Based on Inferred Parameters:
%>                  - Utilize inferred parameters to run forward models.
%>                  - Compare modeled data with experimental data and visualize the results.
%======================================================================================================================
% Execution of Bayesian Analysis
tic
myBayesianAnalysis_Step1 = uq_createAnalysis(BayesOpts_Step1); % Create Bayesian analysis
uq_time = toc; % Record the time taken for the analysis

% Filtering MCMC chains with low acceptance rates
uq_display( myBayesianAnalysis_Step1, 'acceptance', true)
acceptance =  myBayesianAnalysis_Step1.Results.Acceptance;
[~,tolL,tolU,tolC] = isoutlier(acceptance, 'ThresholdFactor', 2);
TF=acceptance<min(max(tolL,0.1),tolC);
badchains = find(TF);
yline(tolL,'b--');
yline(tolU,'b--');
yline(tolC,'r--');
uq_postProcessInversion(myBayesianAnalysis_Step1,'badChains',badchains,'pointEstimate','MAP');
hold on
scatter(badchains, acceptance(badchains), 'red', 'filled')
legend('All chains', 'lb', 'ub', 'median')

% Reporting and Displaying Bayesian Analysis Results
uq_print(myBayesianAnalysis_Step1); % Print out the analysis results
uq_display(myBayesianAnalysis_Step1); % Display results (prior and posterior distributions)

for j =1:10
    fig = figure(j);
    filename = sprintf(strcat('SaveResults_BPM_Step1_Figure%d.jpg'),j);
    saveas(fig,filename,'jpg');
    close(figure(j))
end

% Extracting optimal parameters from the Bayesian analysis
optpar_Step1_log = myBayesianAnalysis_Step1.Results.PostProc.PointEstimate.X{1, 1};
optpar_Step1 = 10.^optpar_Step1_log(1:6); % Transforming the parameters back from log scale

% Saving the Bayesian analysis results
filename_mat_Step1 = sprintf('SaveResults_BPM_Step1_%s.mat',ModelName);
save(filename_mat_Step1); % Save the current results

% Forward Modeling using Inferred Parameters
for i = 1:NumCurve
    Dm = optpar_Step1(3);
    Sc = mu/(rho_g*Dm);
    Re = v*dp*rho_g/mu;
    Dz = Dm * (20+ 0.5*Sc*Re)/epsilon;
    kg = Dm/dp*(2.0+ 1.8 * Re^0.5 * Sc^(1/3));
    Rp = dp/2;
    De = optpar_Step1(4);
    
    if ModelNum==1
        q0star =Isotherm_Langmuir(cFeed(i),optpar_Step1(1),optpar_Step1(2));
    end
    kc= optpar_Step1(5);
    kd= optpar_Step1(6);
    inv_K = (Rp*rho_b*q0star)./(3*kg*cFeed(i)*epsilon) + (Rp^2*rho_b*q0star)./(15*De*cFeed(i)*epsilon);
    K_G = 1/inv_K;
    par = [K_G,optpar_Step1(1),optpar_Step1(2),Dz,optpar_Step1(5),optpar_Step1(6)]; %disp(par)
    
    DisplayMAP_Step1= sprintf('Optimal Parameter Set of %d ppm (Step 1)',cFeed_ppm(i));
    disp(DisplayMAP_Step1);
    Optpar_Save_Step1(i,:) = [optpar_Step1 K_G Dz];
    disp(Optpar_Save_Step1(i,:));
    
    filename = sprintf('Ortiz_error_%d.csv',cFeed_ppm(i));
    exp_inf= csvread(filename,1,0);
    t_exp_plot = exp_inf(:,1);
    c_exp_plot = exp_inf(:,2);
    tf_plot = 60*max(t_exp_plot);
    [t,c,~]= Model_DualFunc_GasPurification(dt,ngrid,dz,par,epsilon,cFeed(i),tf_plot,ModelNum,v,rho);
    c_breakthrough=(c(ngrid,:)+c(ngrid+1,:))./2;
    c_est =  c_breakthrough/convert;
    c_est = c_est/(1e-3 /34.08);
    
    fig=figure(1)
    hold on
    plot(t_exp_plot,c_exp_plot,'o')
    hold on
    plot(t/60,c_est)
    xlabel('t [min]')
    ylabel('c [ppm]')
end
filename = sprintf('SaveResults_BPM_Step1_BreakthroughCurve.jpg');
saveas(fig,filename,'jpg');
close(figure(1))

% Saving the results of forward modeling
filename_csv_Step1 = sprintf('BPM_Step1_OptimalParameter_%s.csv',ModelName);
header = {'KL', 'qm','Dm', 'De','kc','kd','k','Dz'};
fid = fopen(filename_csv_Step1, 'w');
fprintf(fid, '%s,', header{1:end-1});
fprintf(fid, '%s\n', header{end});
fclose(fid);
dlmwrite(filename_csv_Step1, Optpar_Save_Step1, '-append', 'precision', '%.11f', 'delimiter', ','); % Append results to a CSV file
%% Step 2 of Two-Step Bayesian Framework
%% Initial Preparation for Parameter Estimation
%==========================================================================================================================================================================================================================================
%> @details     This section of BPM_BPM_DualFunc_GasPurification_Main.m encompasses the preparation of experimental data and the initial setup for Bayesian analysis in dual-functional gas purification models. The process includes:

%>              (1) Data Analysis and Preparation for Bayesian Inference:
%>                  - Performing linear interpolation, numerical integration, and normalization of data.
%>                  - Preparing data for Bayesian inference by interpolating and smoothing concentration data.
%>                  - Compiling differential concentration rate data and normalized time-to-breakthrough data.

%>              (2) Bayesian Analysis Setup:
%>                  - Defining the forward model for the Bayesian analysis using a custom function.
%>                  - Setting prior distributions for various parameters based on expected ranges.
%>                  - Configuring the Bayesian inversion solver, including the selection of a sampler and defining the number of chains and steps.
%>                  - Preparing the discrepancy model and options for the Bayesian analysis.

%>              This comprehensive approach ensures that all the necessary steps are taken to accurately prepare the experimental data and set up the Bayesian analysis for effective parameter identification in dual-functional gas purification models
%===========================================================================================================================================================================================================================================================
disp('Step 2 of 2-step BPM for Dual-Functional Gas Purification')
index = 1;

post_samples_step1=myBayesianAnalysis_Step1.Results.PostProc.PostSample(:,1:4,:);

for j =1:4
    post_step1(:,j ) = reshape(post_samples_step1(:,j,:), [], 1);
end

for j = 1:4
    data_step1 = post_step1(:, j); % log-scale 데이터 열 추출
    data_step1 = data_step1(:);
    
    Q1 = prctile(data_step1, 25);
    Q3 = prctile(data_step1, 75);
    IQR = Q3 - Q1;
    
    Low(j)  = Q1 - 1.5 * IQR;
    Upp(j) = Q3 + 1.5 * IQR;
    
    step1_bounds_X(1, j) = Low(j);
    step1_bounds_X(2, j) = Upp(j);
end

% Processing each experimental data curve
filename = sprintf('Ortiz_error_%d.csv',cFeed_ppm(index));
exp_inf= csvread(filename,1,0);

t_exp_tmp = exp_inf(:,1);
c_exp_tmp = exp_inf(:,2);

idx_t= find(t_exp_tmp~=0);
t_exp_tmp = t_exp_tmp(idx_t);
c_exp_tmp = c_exp_tmp(idx_t);

idx_tmax=find(t_exp_tmp == max(t_exp_tmp),1,'first');

t_data_tmp = t_exp_tmp(1:idx_tmax)*60; %[s]
c_convert = exp_inf(1:idx_tmax,2);
c_convert2 = c_convert * convert;
c_convert2 = c_convert2 * 1e-3 /34.08;  %[mol/m^3]

t_data = linspace(t_data_tmp(1),max(t_data_tmp),3*Ndata);% Time vector
c_data = interp1(t_data_tmp,c_convert2/cFeed(index),t_data,'linear');

idx1=find(c_data <= 0,1,'last');
dt_tmp = t_data(2:end)-t_data(1:end-1);
[max_dcdt_tmp,idx3] = max(diff(c_data(1:end))./dt_tmp);

t_interp = linspace(t_data(idx1),t_data(end),Ndata);% Time vector
c_interp = interp1(t_data(idx1:end),c_data(idx1:end)*cFeed(index),t_interp, 'linear');
tf(index)= t_interp(end);

hold on
plot(t_interp/60,c_interp/cFeed(index), 'ro', 'LineWidth', 1.5)

% Finalizing the output data for Bayesian analysis : Normalize data
max_c = max(c_interp);
data_c = c_interp/max_c;
dcdt_tmp = data_c;
t_exp = t_interp;
c_exp = c_interp;
myData.y = dcdt_tmp; disp(size(myData.y))

% Setting up the forward model for Bayesian analysis
ModelOpts_Step2.mHandle = @(par) uq_DualFunc_GasPurification_Step2(par,t_exp,c_exp,dt,ngrid,dz,cFeed,epsilon,tf,v,ModelNum,rho,rho_b,rho_g,dp,mu,max_c,index);
ModelOpts_Step2.isVectorized = true;
myForwardModel_Step2 = uq_createModel(ModelOpts_Step2);
BayesOpts_Step2.ForwardModel = myForwardModel_Step2;

% Defining prior distributions for Bayesian parameters
PriorOpts_Step2.Marginals(1).Name = 'b1';
PriorOpts_Step2.Marginals(1).Type = 'Uniform';
PriorOpts_Step2.Marginals(1).Parameters = step1_bounds_X(:,1);

PriorOpts_Step2.Marginals(2).Name = 'qm';
PriorOpts_Step2.Marginals(2).Type = 'Uniform';
PriorOpts_Step2.Marginals(2).Parameters = step1_bounds_X(:,2);

PriorOpts_Step2.Marginals(3).Name = 'Dm';
PriorOpts_Step2.Marginals(3).Type = 'Uniform';
PriorOpts_Step2.Marginals(3).Parameters = step1_bounds_X(:,3);

PriorOpts_Step2.Marginals(4).Name = 'De';
PriorOpts_Step2.Marginals(4).Type = 'Uniform';
PriorOpts_Step2.Marginals(4).Parameters = step1_bounds_X(:,4);

PriorOpts_Step2.Marginals(5).Name = 'kc';
PriorOpts_Step2.Marginals(5).Type = 'Uniform';
PriorOpts_Step2.Marginals(5).Parameters = [-6 -3];

PriorOpts_Step2.Marginals(6).Name = 'kd';
PriorOpts_Step2.Marginals(6).Type = 'Uniform';
PriorOpts_Step2.Marginals(6).Parameters = [-6 -3];
% Creating the prior distribution as an input object
myPriorDist_Step2 = uq_createInput(PriorOpts_Step2);
BayesOpts_Step2.Prior = myPriorDist_Step2;

% Setting up the discrepancy model for the Bayesian analysis
SigmaOpts_Step2.Marginals.Name = 'sigma2';
SigmaOpts_Step2.Marginals.Type = 'Uniform';
SigmaOpts_Step2.Marginals.Parameters = [0 (0.01)^2];
mySigmaDist = uq_createInput(SigmaOpts_Step2);
DiscrepancyOptsUnknownDisc_Step2.Type = 'Gaussian';
DiscrepancyOptsUnknownDisc_Step2.Prior = mySigmaDist;
BayesOpts_Step2.Discrepancy = DiscrepancyOptsUnknownDisc_Step2;

% Configuring the Bayesian inversion solver and sampler
Solver_Step2.Type = 'MCMC';
Solver_Step2.MCMC.Sampler = 'AIES';
Solver_Step2.MCMC.NChains = 100;
Solver_Step2.MCMC.Visualize.Interval = 10;
Solver_Step2.MCMC.Steps =500;
Solver_Step2.MCMC.Visualize.Parameters = [1 2 3 4 5 6];

% Finalizing Bayesian inversion setup
BayesOpts_Step2.Type = 'Inversion';
BayesOpts_Step2.Data = myData;
BayesOpts_Step2.Solver = Solver_Step2;
%% Execution of Bayesian Analysis and Forward Modeling Using IUQ Results
%======================================================================================================================
%> @details     This section performs the Bayesian analysis for gas adsorption parameter identification and conducts forward modeling using the inferred parameters.
%>              The process includes:

%>              (1) Bayesian Analysis Execution:
%>                  - Perform Bayesian analysis using the specified options.
%>
%>              (2) Post-Processing of Bayesian Results:
%>                  - Filter out MCMC chains with low acceptance rates.
%>                  - Display and analyze the acceptance rates of the chains.
%>
%>              (3) Results Reporting and Visualization:
%>                  - Print a detailed report of the Bayesian analysis results.
%>                  - Display results including prior and posterior distributions.
%>
%>              (4) Forward Modeling Based on Inferred Parameters:
%>                  - Utilize inferred parameters to run forward models.
%>                  - Compare modeled data with experimental data and visualize the results.
%======================================================================================================================
% Execution of Bayesian Analysis
tic
myBayesianAnalysis_Step2 = uq_createAnalysis(BayesOpts_Step2); % Create Bayesian analysis
uq_time = toc; % Record the time taken for the analysis

% Filtering MCMC chains with low acceptance rates
uq_display( myBayesianAnalysis_Step2, 'acceptance', true)
acceptance =  myBayesianAnalysis_Step2.Results.Acceptance;
[~,tolL,tolU,tolC] = isoutlier(acceptance, 'ThresholdFactor', 2);
TF=acceptance<min(max(tolL,0.1),tolC);
badchains = find(TF);
yline(tolL,'b--');
yline(tolU,'b--');
yline(tolC,'r--');
uq_postProcessInversion(myBayesianAnalysis_Step2,'badChains',badchains,'pointEstimate','MAP');
hold on
scatter(badchains, acceptance(badchains), 'red', 'filled')
legend('All chains', 'lb', 'ub', 'median')

% Reporting and Displaying Bayesian Analysis Results
uq_print(myBayesianAnalysis_Step2); % Print out the analysis results
uq_display(myBayesianAnalysis_Step2); % Display results (prior and posterior distributions)

for j =1:10
    fig = figure(j);
    filename = sprintf(strcat('SaveResults_BPM_Step1_%dth_curve_Figure%d.jpg'),index,j);
    saveas(fig,filename,'jpg');
    close(figure(j))
end

% Extracting optimal parameters from the Bayesian analysis
optpar_Step2_log = myBayesianAnalysis_Step2.Results.PostProc.PointEstimate.X{1, 1};
optpar_Step2 = 10.^optpar_Step2_log(1:6); % Transforming the parameters back from log scale

% Saving the Bayesian analysis results
filename_mat_Step2 = sprintf('SaveResults_BPM_Step2_%s.mat',ModelName);
save(filename_mat_Step2); % Save the current results

% Forward Modeling using Inferred Parameters
Dm = optpar_Step2 (3);
Sc = mu/(rho_g*Dm);
Re = v*dp*rho_g/mu;
Dz = Dm * (20+ 0.5*Sc*Re)/epsilon;
kg = Dm/dp*(2.0+ 1.8 * Re^0.5 * Sc^(1/3));
Rp = dp/2;
De = optpar_Step2(4);

if num==1
    q0star =Isotherm_Langmuir(cFeed(index),optpar_Step2 (1),optpar_Step2 (2));
end

kc= optpar_Step2(5);
kd= optpar_Step2(6);

inv_K = (Rp*rho_b*q0star)./(3*kg*cFeed(index)*epsilon) + (Rp^2*rho_b*q0star)./(15*De*cFeed(index)*epsilon);
K_G = 1/inv_K;
Isotherm_param = [optpar_tmp(1),optpar_tmp(2)];
par = [K_G,optpar_Step2 (1),optpar_Step2 (2),Dz,kc,kd];

DisplayMAP_Step2= sprintf('Optimal Parameter Set of %d ppm (Step 2)',cFeed_ppm(index));
disp(DisplayMAP_Step2);
Optpar_Save_Step2 = [optpar_Step2 K_G Dz];
disp(Optpar_Save_Step2);

filename = sprintf('Ortiz_error_%d.csv',cFeed_ppm(index));
exp_inf= csvread(filename,1,0);
t_exp_plot = exp_inf(:,1);
c_exp_plot = exp_inf(:,2);
tf_plot = 60*max(t_exp_plot);

[t,c,~]= Model_DualFunc_GasPurification(dt,ngrid,dz,par,epsilon,cFeed(index),tf_plot,num,v,rho);
c_breakthrough=(c(ngrid,:)+c(ngrid+1,:))./2;
c_est = c_breakthrough/convert;
c_est = c_est/(1e-3 /34.08);

fig=figure(1);
hold on
plot(t_exp_plot,c_exp_plot,'o')
hold on
plot(t/60,c_est)
xlabel('t [min]')
ylabel('c [ppm]')
filename = sprintf('SaveResults_BPM_Step2_BreakthroughCurve.jpg');
saveas(fig,filename,'jpg');
close(figure(1))

% Saving the results of forward modeling
filename_csv_Step2 = sprintf('BPM_Step2_OptimalParameter_%s.csv',ModelName);
header = { 'KL', 'qm','Dm', 'De','kc','kd','k','Dz'};
fid = fopen(filename_csv_Step2, 'w');
fprintf(fid, '%s,', header{1:end-1});
fprintf(fid, '%s\n', header{end});
fclose(fid);
dlmwrite(filename_csv_Step2, Optpar_Save_Step2, '-append', 'precision', '%.11f', 'delimiter', ','); % Append results to a CSV file