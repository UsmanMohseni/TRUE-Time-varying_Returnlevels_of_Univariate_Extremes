%-------------------------------------------------------------------------
%                           RUN TRUE (GUI)
%-------------------------------------------------------------------------

%Reference Publication: Mohseni U, Gailakoti M, Vinnarasi R, 2025
% Non-Stationary Multi-Covariate Flood Frequency Analysis with TRUE: 
% A Tool for Time-varying Return levels of Univariate Extremes
% 

%TRUE: Time-varying Return level for Univariate Extremes

% TRUE MATLAB code is developed by Usman Mohseni, Research Scholar, 
% Department of Civil Engineering,
% Indian Institute of Technology Roorkee, Uttarakhand, India.

% For questions and permissions, please contact 
% Dr. Vinnarasi Rajendran: vinnarasi@ce.iitr.ac.in or
% Usman Mohseni: mohseni_ua@ce.iitr.ac.in 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disclaimer: The TRUE (Time-varying Return level for Univariate Extremes) 
% toolbox is made available freely for research and educational purposes. 
% It is distributed without any warranties or guarantees, either expressed 
% or implied, regarding the accuracy, completeness, or suitability of the 
% code, results, or visualizations. Users are fully responsible for how 
% they apply TRUE and for any consequences that may arise from its use. 
% The developers and affiliated institutions shall not be held liable for 
% any loss, damage, or costs (whether direct, indirect, incidental, or 
% consequential) connected to the use of this software. The algorithms and 
% features may be updated or modified without prior notice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc
clear
close all


fileLocation = "E:\TRUE"; % Edit location to store the figures
finalsave = "E:\TRUE\All_Results.mat"; % Edit workspace save location

main();

cd('TRUE_Code_Files');

TRUE_GUI_NSMCFFA

