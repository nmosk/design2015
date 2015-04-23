clear;
clc;

% LEVEL 4: DISTILLATION COLUMN 
% ------------------------------


% MOLECULAR WEIGHTS [g/mol] -----------------------------
MW_eB = 106.17;
MW_St = 104.18;
MW_B = 78.11;
MW_H2 = 2.016;
MW_eN = 28.05;
MW_T = 92.138;
MW_Me = 16.04;
        MW = [106.17 104.18 78.1118 92.1384 28.0532 2.016 16.0425 18.016]; %[g/mol]
        
% DENSITIES [g/mL] -----------------------------
rho_eB = 0.867;
rho_St = 0.906;
rho_B = 0.874;
rho_H2 = 0.0000899;
rho_eN = 1.153;
rho_T = 0.865;
rho_Me = .0006672;
rho_H20 = 1;
        rho = [rho_eB rho_St rho_B rho_T rho_eN rho_H2 rho_Me rho_H20]; %[g/mL]
% -------------------------------------------------------
% -------------------------------------------------------

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


% DESCRIPTION ===========================================
%   Calculates Level 4: Distillation column design
% =======================================================

%% First column B T E | S
disp('Column 1 BTE_S split:')

%%
% SPECIFY THESE TERMS
% --------------------
F = 204; % feed molar flowrate [mol/hr]
q = 1; % fraction of feed that is liquid
zF = [0.25 0.56 0.07 0.12 0 0 0 0]; % composition across all phases
T = 303; % temperature in K
P = 1; % pressure in bar

% the species in the distillate and bottoms
% put a 1 if the species is in, put a 0 if it is not
species_D = [1 0 1 1 0 0 0 0];
species_B = [0 1 0 0 0 0 0 0];

% HK and LK order
% [A B C D]
% indices of the species so that HK_LK(1) = A, HK_LK(2) = B, etc
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene
    
 % *************************************************************************   
 % NOTE: PLEASE SEE rmin CALCULATION SECTION TO CHANGE RELATIVE VOLATILITIES
 %       IF DOING A DIFFERENT SPLIT
 % *************************************************************************
 
    HK_LK = [3 4 1 2];

% CALCULATES D, B, xB, xD
% --------------------
F_D = F.*zF.*species_D;
F_B = F.*zF.*species_B;
D = sum(F_D)
B = sum(F_B)
x_D = F_D./D;
x_B = F_B./B;


% Checks that D+B = F
if D+B > 5+F | D+B < F-5
    disp('D+B does not equal F')
else
    disp('D+B equals F -- CHECK 1')
end

% -------------------- 
%%
% CALCULATE MINIMUM REFLUX AND REFLUX
% use Doherty's book, chapter 4 for pseudo 4-composition mixture if >4
% species
% --------------------

% calculating relative volatilities of multi-component mixture
% assuming that volatility is constant
% alpha = (y_LK/x_LK) / (y_HK/x_HK) = K_LK/K_HK

    [RelVol] = RelVol_func(T); % where St is the reference component
    % for an AB/CD split, A = LK and C = HK
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene

 RelVol(1) = RelVol(2)*RelVol(3)*RelVol(1)
 RelVol(2) = RelVol(2)*RelVol(3)
    
% minimum reflux for an AB/CD split
% r_min = ( (RelVol(3)*zF(HK_LK(1))/ (RelVol(1)-RelVol(3))) + ((RelVol(3)*(zF(HK_LK(2))+zF(HK_LK(3))))/(RelVol(2)-RelVol(3))) ) / ((zF(HK_LK(1))+zF(HK_LK(2))) * (1+(zF(HK_LK(1))*(zF(HK_LK(3))+zF(HK_LK(4))))));

% minimum reflux for an ABC/D split
r_min =  ((zF(HK_LK(1))/(RelVol(1)-1)) + (zF(HK_LK(2))/(RelVol(2)-1)) +((zF(HK_LK(3))+zF(HK_LK((4)))/(RelVol(3)-1)))) / ((1-zF(HK_LK((4))))*(1+(zF(HK_LK(4))*(zF(HK_LK(1))+zF(HK_LK(2))))));

    % ------------

r = r_min*1.5
% --------------------

% *** after this point, no more choice in reboiler ratio

s = (D/B)*(r+q)-(1-q); % (eqn 3.35)
%%
% CALCULATE MINIMUM NUMBER OF STAGES
% use Fenske equation(eqn 4.16)
fHK_B = 0.995;
fLK_D = 0.995;
N_min_small = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(3));

fHK_B = 0.999;
fLK_D = 0.999;
N_min_big = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(3));

N_min = N_min_small;
% --------------------

% CALCULATE THEORETICAL AND REAL NUMBER OF STAGES
% use FUG method (p. 136)
% --------------------

% calculating N theoretical

% Function that calculates N_theory values for r_min to 1
% uses the Gilliland design method equation (eqn 4.56)
% also outputs a plot of (N-N_min)/(N+1) vs (r-r_min)/(r+1)
% [ N_theory_array, r_array, N_theory ] = Ntheory_func( r_min, N_min, r );
[ N_theory ] = Ntheory_func( r_min, N_min, r );
N_real = 2.*double(N_theory)

% % O'CONNELL CORRELATION p. 260 (eqn 6.2)
% a = 0.24; 
% mu % viscosity of the liquid mixture at the feed composition evaluated at Tavg and Pavg in column
% alpha % volatility between the key components evaluated at Tavg and Pavg in column
% % The relative volatility is determined for the 2 key components at average column conditions
% mu_0 = 10^-3 % [Pa*s] aka (1 centipoise)
% N_real = N_theory./(exp(-sqrt(alpha*mu/mu_0))*(1-a_param)+a)^-1 % (eqn 6.2)
% --------------------

%%
% CALCULATE VAPOR RATES [mol/hr] IN COLUMN
% --------------------
% ****************> these need to equal each other if q = 1 !!!!!!!!!
v_B = s*B; % in bottoms
v_T = (r+1)*D; % in tops

% CROSS-CHECK: 

if v_B-v_T < 0.001 && v_B-v_T > -.001  %(eqn 3.39)
disp('cross check of vB = vT passed!')
else
    disp('cross check failed')
end
    
%              when q = 1, v_B = v_T = V
V=v_B
% --------------------

%%
% CALCULATE COST MODEL
% shortcut for comparison to different split methods 
% --------------------
[ cost_BTE_S_split, cost_check1 ] = CostModel_func(N_real, V)
% --------------------
%% Begin second column B/TE
disp('Column 2, B/TE:')

%%
% SPECIFY THESE TERMS
% --------------------
F = F; % feed molar flowrate [mol/hr]
q = 1; % fraction of feed that is liquid
zF = [0.25 0.56 0.07 0.12 0 0 0 0]; % composition across all phases
T = 303; % temperature in K
P = 1; % pressure in bar

% the species in the distillate and bottoms
% put a 1 if the species is in, put a 0 if it is not
species_D = [0 0 1 0 0 0 0 0];
species_B = [1 0 0 1 0 0 0 0];

% HK and LK order
% [A B C D]
% indices of the species so that HK_LK(1) = A, HK_LK(2) = B, etc
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene
    
 % *************************************************************************   
 % NOTE: PLEASE SEE rmin CALCULATION SECTION TO CHANGE RELATIVE VOLATILITIES
 %       IF DOING A DIFFERENT SPLIT
 % *************************************************************************
 
    HK_LK = [3 4 1 2];

% CALCULATES D, B, xB, xD
% --------------------
F_D = F.*zF.*species_D;
F_B = F.*zF.*species_B;
D = sum(F_D)
B = sum(F_B)
x_D = F_D./D;
x_B = F_B./B;
new_F = D + B;
xF = (F_D + F_B)./new_F;

% Checks that D+B = F
if D+B > 5+new_F | D+B < new_F-5
    disp('D+B does not equal F')
else
    disp('D+B equals F -- CHECK 1')
end
% --------------------
%%
% CALCULATE MINIMUM REFLUX AND REFLUX
% use Doherty's book, chapter 4 for pseudo 4-composition mixture if >4
% species
% --------------------

% calculating relative volatilities of multi-component mixture
% assuming that volatility is constant
% alpha = (y_LK/x_LK) / (y_HK/x_HK) = K_LK/K_HK

    [RelVol] = RelVol_func(T); % where St is the reference component
    % for an AB/CD split, A = LK and C = HK
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene
 %   RelVol = [# # # 1]; % relative volatilities [A B C D] 
    % where D=1 since St is the reference component

% minimum reflux for an AB/CD split
% r_min = ( (RelVol(3)*zF(HK_LK(1))/ (RelVol(1)-RelVol(3))) + ((RelVol(3)*(zF(HK_LK(2))+zF(HK_LK(3))))/(RelVol(2)-RelVol(3))) ) / ((zF(HK_LK(1))+zF(HK_LK(2))) * (1+(zF(HK_LK(1))*(zF(HK_LK(3))+zF(HK_LK(4))))));

% minimum reflux for an ABC/D split
% r_min =  ((zF(HK_LK(1))*RelVol(1)) + (zF(HK_LK(2))*RelVol(2)) +((zF(HK_LK(3))+zF(HK_LK((4)))*RelVol(3)))) / ((1-zF(HK_LK((4))))*(1+(zF(HK_LK(4))*(zF(HK_LK(1))+zF(HK_LK(2))))));

% minimum reflux for an A/BCD split
% r_min = ((RelVol(2)*(zF(HK_LK(1))+zF(HK_LK(2))))/(zF(HK_LK(1))*(RelVol(1)-RelVol(2)))) + (RelVol(3)*zF(HK_LK(3))/(zF(HK_LK(1))*(RelVol(1)-RelVol(3)))) + (zF(HK_LK(4))/(RelVol(1)-1));

% minimum reflux for an A/BC split
 r_min = ( (RelVol(2)*(xF(HK_LK(1))+xF(HK_LK(2)))) / (xF(HK_LK(1))*(RelVol(1)*RelVol(2)-RelVol(2))) ) + ( xF(HK_LK(3)) / (xF(HK_LK(1))*(RelVol(1)*RelVol(2)-1)) );

% minimum reflux for an AB/C split
% r_min = ( ((zF(HK_LK(2))+zF(HK_LK(3)))/(RelVol(2)-1)) + (zF(HK_LK(1))/(RelVol(1)-1)) ) / ( (zF(HK_LK(1))+zF(HK_LK(2)))*(1+(zF(HK_LK(1))*zF(HK_LK(3)))) );

% FOR BINARY -
%     % alpha = (y_i/x_i) / (y_j/x_j) = K_i/K_j
%     alpha = (y(1)/x(1)) / (y(2)/x(2)); 
%     % where i is the more volatile component and j is the less volatile
%     % component
%     r_min = (alpha(1)-1)^-1 * ( (x_D/z_F)-(alpha(1)*(1-x_D)/(1-z_F)) ); %
%     % Underwood equation for a binary mixture from "Distillation: Fundamentals
%     % and Principles" edited by Andrzej Gorak and Eva Sorensen (eqn 4.71)
    % ------------

r = r_min*1.5
% --------------------

% *** after this point, no more choice in reboiler ratio

s = (D/B)*(r+q)-(1-q); % (eqn 3.35)

%%
%%
% CALCULATE MINIMUM NUMBER OF STAGES
% use Fenske equation(eqn 4.16)
fHK_B = 0.995;
fLK_D = 0.995;
N_min_small = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(1));

fHK_B = 0.999;
fLK_D = 0.999;
N_min_big = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(1));

N_min = N_min_small;
% --------------------

% CALCULATE THEORETICAL AND REAL NUMBER OF STAGES
% use FUG method (p. 136)
% --------------------

% calculating N theoretical

% Function that calculates N_theory values for r_min to 1
% uses the Gilliland design method equation (eqn 4.56)
% also outputs a plot of (N-N_min)/(N+1) vs (r-r_min)/(r+1)
% [ N_theory_array, r_array, N_theory ] = Ntheory_func( r_min, N_min, r );
[ N_theory ] = Ntheory_func( r_min, N_min, r );
N_real = 2.*double(N_theory)

% % O'CONNELL CORRELATION p. 260 (eqn 6.2)
% a = 0.24; 
% mu % viscosity of the liquid mixture at the feed composition evaluated at Tavg and Pavg in column
% alpha % volatility between the key components evaluated at Tavg and Pavg in column
% % The relative volatility is determined for the 2 key components at average column conditions
% mu_0 = 10^-3 % [Pa*s] aka (1 centipoise)
% N_real = N_theory./(exp(-sqrt(alpha*mu/mu_0))*(1-a_param)+a)^-1 % (eqn 6.2)
% --------------------


%%
% CALCULATE VAPOR RATES [mol/hr] IN COLUMN
% --------------------
% ****************> these need to equal each other if q = 1 !!!!!!!!!
v_B = s*B; % in bottoms
v_T = (r+1)*D; % in tops

% CROSS-CHECK: 

if v_B-v_T < 0.001 && v_B-v_T > -.001  %(eqn 3.39)
disp('cross check of vB = vT passed!')
else
    disp('cross check failed')
end
    
%              when q = 1, v_B = v_T = V
V=v_B
% --------------------

%%
% CALCULATE COST MODEL
% shortcut for comparison to different split methods 
% --------------------
[ cost_B_TE_split, cost_check2 ] = CostModel_func(N_real, V)
% --------------------
%% Begin column 3 T/E
disp('Column 3 T/E:')

%%
% SPECIFY THESE TERMS
% --------------------
F = 204; % feed molar flowrate [mol/hr]
q = 1; % fraction of feed that is liquid
zF = [0.25 0.56 0.07 0.12 0 0 0 0]; % composition across all phases
T = 303; % temperature in K
P = 1; % pressure in bar

% the species in the distillate and bottoms
% put a 1 if the species is in, put a 0 if it is not
species_D = [0 0 0 1 0 0 0 0];
species_B = [1 0 0 0 0 0 0 0];

% HK and LK order
% [A B C D]
% indices of the species so that HK_LK(1) = A, HK_LK(2) = B, etc
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene
    
 % *************************************************************************   
 % NOTE: PLEASE SEE rmin CALCULATION SECTION TO CHANGE RELATIVE VOLATILITIES
 %       IF DOING A DIFFERENT SPLIT
 % *************************************************************************
 
    HK_LK = [3 4 1 2];

% CALCULATES D, B, xB, xD
% --------------------
F_D = F.*zF.*species_D;
F_B = F.*zF.*species_B;
D = sum(F_D)
B = sum(F_B)
x_D = F_D./D;
x_B = F_B./B;
new_F = D + B;
xF = (F_D + F_B)./new_F;

% Checks that D+B = F
if D+B > 5+new_F | D+B < new_F-5
    disp('D+B does not equal F')
else
    disp('D+B equals F -- CHECK 1')
end

% --------------------
%%
% CALCULATE MINIMUM REFLUX AND REFLUX
% use Doherty's book, chapter 4 for pseudo 4-composition mixture if >4
% species
% --------------------

% calculating relative volatilities of multi-component mixture
% assuming that volatility is constant
% alpha = (y_LK/x_LK) / (y_HK/x_HK) = K_LK/K_HK

    [RelVol] = RelVol_func(T); % where St is the reference component
    % for an AB/CD split, A = LK and C = HK
    % for an ABC/D split, C = LK and D = HK
    % eB is the LK and St is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Benzene , B=Toluene, C=Ethylbenzene, D=Styrene
 %   RelVol = [# # # 1]; % relative volatilities [A B C D] 
    % where D=1 since St is the reference component

% minimum reflux for an AB/CD split
% r_min = ( (RelVol(3)*zF(HK_LK(1))/ (RelVol(1)-RelVol(3))) + ((RelVol(3)*(zF(HK_LK(2))+zF(HK_LK(3))))/(RelVol(2)-RelVol(3))) ) / ((zF(HK_LK(1))+zF(HK_LK(2))) * (1+(zF(HK_LK(1))*(zF(HK_LK(3))+zF(HK_LK(4))))));

% minimum reflux for an ABC/D split
% r_min =  ((zF(HK_LK(1))*RelVol(1)) + (zF(HK_LK(2))*RelVol(2)) +((zF(HK_LK(3))+zF(HK_LK((4)))*RelVol(3)))) / ((1-zF(HK_LK((4))))*(1+(zF(HK_LK(4))*(zF(HK_LK(1))+zF(HK_LK(2))))));

% minimum reflux for an A/BCD split
% r_min = ((RelVol(2)*(zF(HK_LK(1))+zF(HK_LK(2))))/(zF(HK_LK(1))*(RelVol(1)-RelVol(2)))) + (RelVol(3)*zF(HK_LK(3))/(zF(HK_LK(1))*(RelVol(1)-RelVol(3)))) + (zF(HK_LK(4))/(RelVol(1)-1));

% minimum reflux for an A/BC split
% r_min = ( (RelVol(2)*(zF(HK_LK(1))+zF(HK_LK(2)))) / (zF(HK_LK(1))*(RelVol(1)-RelVol(3))) ) + ( zF(HK_LK(3)) / (zF(HK_LK(1))*(RelVol(1)-1)) );

% minimum reflux for an AB/C split
% r_min = ( ((zF(HK_LK(2))+zF(HK_LK(3)))/(RelVol(2)-1)) + (zF(HK_LK(1))/(RelVol(1)-1)) ) / ( (zF(HK_LK(1))+zF(HK_LK(2)))*(1+(zF(HK_LK(1))*zF(HK_LK(3)))) );

%FOR BINARY -
    % where i is the more volatile component and j is the less volatile
    % component
    r_min = (RelVol(2)-1)^-1 * ((x_D(HK_LK(2))/xF(HK_LK(2))- RelVol(2)*(1-x_D(HK_LK(2))/(1-xF(HK_LK(2)))))); % (3.57) eqn Doherty

r = r_min*1.5
% --------------------

% *** after this point, no more choice in reboiler ratio

s = (D/B)*(r+q)-(1-q); % (eqn 3.35)

%%
% CALCULATE MINIMUM NUMBER OF STAGES
% use Fenske equation(eqn 4.16)
fHK_B = 0.995;
fLK_D = 0.995;
N_min_small = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(2));

fHK_B = 0.999;
fLK_D = 0.999;
N_min_big = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(2));

N_min = N_min_small;
% --------------------

% CALCULATE THEORETICAL AND REAL NUMBER OF STAGES
% use FUG method (p. 136)
% --------------------

% calculating N theoretical

% Function that calculates N_theory values for r_min to 1
% uses the Gilliland design method equation (eqn 4.56)
% also outputs a plot of (N-N_min)/(N+1) vs (r-r_min)/(r+1)
% [ N_theory_array, r_array, N_theory ] = Ntheory_func( r_min, N_min, r );
[ N_theory ] = Ntheory_func( r_min, N_min, r );
N_real = 2.*double(N_theory)

% % O'CONNELL CORRELATION p. 260 (eqn 6.2)
% a = 0.24; 
% mu % viscosity of the liquid mixture at the feed composition evaluated at Tavg and Pavg in column
% alpha % volatility between the key components evaluated at Tavg and Pavg in column
% % The relative volatility is determined for the 2 key components at average column conditions
% mu_0 = 10^-3 % [Pa*s] aka (1 centipoise)
% N_real = N_theory./(exp(-sqrt(alpha*mu/mu_0))*(1-a_param)+a)^-1 % (eqn 6.2)
% --------------------

% CALCULATE VAPOR RATES [mol/hr] IN COLUMN
% --------------------
% ****************> these need to equal each other if q = 1 !!!!!!!!!
v_B = s*B; % in bottoms
v_T = (r+1)*D; % in tops

% CROSS-CHECK: 

if v_B-v_T < 0.001 && v_B-v_T > -.001  %(eqn 3.39)
disp('cross check of vB = vT passed!')
else
    disp('cross check failed')
end
    
%              when q = 1, v_B = v_T = V
V=v_B
% --------------------

%%
% CALCULATE COST MODEL
% shortcut for comparison to different split methods 
% --------------------
[ cost_T_E_split, cost_check3 ] = CostModel_func(N_real, V)
% --------------------

% total_cost = cost_BTE_S_split + cost_B_TE_split + cost_T_E_split
% Cost estimation of all columns
total_cost = cost_check1(1) + cost_check2(1) + cost_check3(1)

%%

