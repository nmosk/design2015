function [ installed_cost_column, purchased_cost_heatex, installed_cost_heatex ] = EquipCosts_func( height_column, diameter_column, A )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% INSTALLED COST OF COLUMN
% Douglas Appendix E
% --------------------
    % Correction factors
    % from (table E-2-6)
    F_s = 2.2; % for 12 inch tray spacing
    F_t = 0; % for plate tray
    F_m = 1.7;  % for stainless steel
    F_c = F_s+F_t+F_m;
    
    He = height_column*3.28; % height of column in [ft]
    Di = diameter_column*3.28; % diameter of column in [ft]
    MAS = 1600; 
    
installed_cost_column = (MAS/280) * 4.7 * Di^1.55*He*F_c;
% --------------------
%%
% PURCHASED AND INSTALLED COST OF HEAT EXCHANGER
% Douglas Appendix E
% --------------------
    % Correction factors
    % from (table E-2-6)
    F_d = 0.85; % for Utube
    F_p = 0.25; % for up to 400 psi
    F_m = 2.81;  % for CS/SS
    F_c = (F_d+F_p)*F_m; 

    A = A*(3.28)^2; % area of heat exchanger in [ft^2]

purchased_cost_heatex = (MAS/280)*(101.3*A^0.65*F_c);
installed_cost_heatex = (MAS/280)*(101.3*A^0.65*(F_c+2.29));
end

