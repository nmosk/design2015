% ABBREVIATIONS of SPECIES ------------------------------
% ethylbenzene -> eB 
% styrene -> St 
% benzene -> B 
% toulene -> T 
% ethylene -> eN 
% hydrogen -> H2 
% methane -> Me 
% water-> H20
% -------------------------------------------------------

% NEW ORDER OF SPECIES IN ARRAYS ----------------------
% eb  -> 1
% St  -> 2
% B  ->  3
% T  ->  4
% eN  -> 5
% H2  -> 6
% Me  -> 7
% H20 -> 8
% -------------------------------------------------------

function [ lambda_D, lambda_B ] = HeatVap_func(x_D,x_B)

% DESCRIPTION ===========================================
%   Calculates latent heat of vaporization of distillate
%   and bottoms
% =======================================================

% ******************************************************
% NOTE: x_D, x_B needs to be arrays of size (1,8).     *
%       if a species does not exist in D or B, its cell*
%       in the array must be '0'                       *
% ******************************************************


% MOLECULAR WEIGHTS [g/mol] -----------------------------
MW_eB = 106.17;
MW_St = 104.18;
MW_B = 78.11;
MW_H2 = 2.016;
MW_eN = 28.05;
MW_T = 92.138;
MW_Me = 16.04;
        MW = [106.17 104.18 78.1118 92.1384 28.0532 2.016 16.0425 18.016]; %[g/mol]
% -------------------------------------------------------

% LATENT HEAT OF VAPORIZATION [J/g] -----------------
% all in [kJ/kg]
lambda_eB = 335;
lambda_St = 363;
lambda_B = 390;
lambda_T = 351; 
lambda_eN = 482;
lambda_H2 = 461;
lambda_Me = 112;
lambda_H20 = 2257;
        lambdas = [lambda_eB lambda_St lambda_B lambda_T lambda_eN lambda_H2 lambda_Me lambda_H20]; % [kJ/kg]
        lambdas = lambdas.*MW; % change to [kJ/kmol]
% -------------------------------------------------------


% FIND LAMBDAS using the arithmatic mean

    % LAMBDA FOR DISTILLATE
    x_D; % composition of distillate
            x_D_elimzeros = x_D; 
            x_D_elimzeros(x_D_elimzeros==0) = [];
            
    lambda_D = sum(x_D.*lambdas)./length(x_D_elimzeros); % [J/mol]

    % LAMBDA FOR BOTTOMS
    x_B; % composition of distillate
            x_B_elimzeros = x_B; 
            x_B_elimzeros(x_B_elimzeros==0) = [];    
            
    lambda_B = sum(x_B.*lambdas)./length(x_B_elimzeros); % [J/mol]

end

