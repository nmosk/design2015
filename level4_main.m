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

%%
% SPECIFY THESE TERMS
% --------------------
% F = 204; % feed molar flowrate [mol/hr]
 F = 90 ; % feed of second distillation column
% F= 77 ; % feed of 3rd distil. 

q = 1; % fraction of feed that is liquid

% zF = [0.25 0.56 0.06 0.12 0 0 0 0]; % composition across all phases
 zF= [.58 0 .14 .27 0 0 0 0]; % composition for 2nd distil. 
% zF= [0.68 0 0 0.31 0 0 0 0] ; % comp for 3rd distil. 

% T = 100 + 273 ; % temperature in K
 T = 81  + 273 ; % temp of 2nd distill column
% T = 131 + 273 ; % temp of 3rd distill column

P = 1; % pressure in bar


% the species in the distillate and bottoms
% put a 1 if the species is in, put a 0 if it is not
% species_D = [1 0 1 1 0 0 0 0]; % column 1
% species_B = [0 1 0 0 0 0 0 0]; % column 1

species_D = [0 0 1 0 0 0 0 0]; % column 2
species_B = [1 0 0 1 0 0 0 0]; % column 2

%  species_D = [0 0 0 1 0 0 0 0]; % column 3
%  species_B = [1 0 0 0 0 0 0 0]; % column 3

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
 
   % HK_LK = [3 4 1 2]; % for column 1
    HK_LK = [3 4 1]; % for column 2
   % HK_LK = [4 1]; % for column 3 

    %%
    
% CALCULATES D, B, xB, xD
% --------------------
F_D = F.*zF.*species_D;
F_B = F.*zF.*species_B;
D = sum(F_D)
B = sum(F_B)
x_D = F_D./D
x_B = F_B./B


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
 %   RelVol = [# # # 1]; % relative volatilities [A B C D] 
    % where D=1 since St is the reference component

% minimum reflux for an AB/CD split
% r_min = ( (RelVol(3)*zF(HK_LK(1))/ (RelVol(1)-RelVol(3))) + ((RelVol(3)*(zF(HK_LK(2))+zF(HK_LK(3))))/(RelVol(2)-RelVol(3))) ) / ((zF(HK_LK(1))+zF(HK_LK(2))) * (1+(zF(HK_LK(1))*(zF(HK_LK(3))+zF(HK_LK(4))))));

% minimum reflux for an ABC/D split **
% r_min =  ((zF(HK_LK(1))*RelVol(1)) + (zF(HK_LK(2))*RelVol(2)) +((zF(HK_LK(3))+zF(HK_LK((4)))*RelVol(3)))) / ((1-zF(HK_LK((4))))*(1+(zF(HK_LK(4))*(zF(HK_LK(1))+zF(HK_LK(2))))));

% minimum reflux for an A/BCD split
% r_min = ((RelVol(2)*(zF(HK_LK(1))+zF(HK_LK(2))))/(zF(HK_LK(1))*(RelVol(1)-RelVol(2)))) + (RelVol(3)*zF(HK_LK(3))/(zF(HK_LK(1))*(RelVol(1)-RelVol(3)))) + (zF(HK_LK(4))/(RelVol(1)-1));

% minimum reflux for an A/BC split ** 
 r_min = ( (RelVol(2)*(zF(HK_LK(1))+zF(HK_LK(2)))) / (zF(HK_LK(1))*(RelVol(1)-RelVol(3))) ) + ( zF(HK_LK(3)) / (zF(HK_LK(1))*(RelVol(1)-1)) );

% minimum reflux for an AB/C split
% r_min = ( ((zF(HK_LK(2))+zF(HK_LK(3)))/(RelVol(2)-1)) + (zF(HK_LK(1))/(RelVol(1)-1)) ) / ( (zF(HK_LK(1))+zF(HK_LK(2)))*(1+(zF(HK_LK(1))*zF(HK_LK(3)))) );

% FOR BINARY -
%     % alpha = (y_i/x_i) / (y_j/x_j) = K_i/K_j
%     alpha = (HK_LK(1)/HK_LK(1)) / (HK_LK(2)/HK_LK(2)); 
%     alpha = 1.1034/1.172
%     % where i is the more volatile component and j is the less volatile
%     % component
%   r_min = (alpha(1)-1)^-1 * ( (x_D/zF)-(alpha(1)*(1-x_D)/(1-zF)) ); %
%     % Underwood equation for a binary mixture from "Distillation: Fundamentals
%     % and Principles" edited by Andrzej Gorak and Eva Sorensen (eqn 4.71)
    % ------------

r = r_min*1.5
% --------------------

% *** after this point, no more choice in reboiler ratio

s = (D/B)*(r+q)-(1-q) % (eqn 3.35)

%%
%%
% CALCULATE MINIMUM NUMBER OF STAGES
% use Fenske equation(eqn 4.16)
fHK_B = 0.995;
fLK_D = 0.995;
N_min_small = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(3))

fHK_B = 0.999;
fLK_D = 0.999;
N_min_big = log((fLK_D/(1-fLK_D))*(fHK_B/(1-fHK_B)))/log(RelVol(3))

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
[ costs, cost_check ] = CostModel_func(N_real, V);
% --------------------



%%
% CALCULATE COST MODEL
% shortcut for comparison to different split methods 
% --------------------
[ costs, cost_check ] = CostModel_func(N_real, V);
% --------------------

%%
% CALCULATE HEAT LOADS (saturated liquid products)
% lambda is the latent heat of vaporization
% estimate by taking an arithmatic average
% --------------------
[ lambda_D, lambda_B ] = HeatVap_func(x_D,x_B);

Q_c = lambda_D*v_T % [J/hr]
Q_r = lambda_B*v_B % [J/hr]
% --------------------

%%

% CALCULATE COLUMN DIMENSIONS
%   height_theory is approximately 30 to 60 cm (12 to 24 inches) but can be
%   smaller or larger for some columns
    
% diameter with (eqn. 6.12)
% -------------------
% CALCULATING HEIGHT
% --------------------------------------------
%             (TABLE 6.1)    TRAY SPACING
%            0.31 m      0.46 m     0.61 m
%             12 in       18 in      24 in
% --------------------------------------------
% c_o (m/h)    252        329         439
%     (ft/h)   828        1080        1440
% c1           2.0         2.3         2.5
% c2           1.0         1.1         1.2
% --------------------------------------------

        height_tray = 0.5;   % spacing between trays [m]
        
        M_v = x_D(1)*MW_eB + x_D(3)*MW_B + x_D(4)*MW_T % molecular weight of the vapor     
        M_l = MW_St % molecular weight of the liquid
        
        
        c_o = 439; % [m/h] % from (table 6.1)
        height_column_min = 3*height_tray;
        height_column = height_column_min+(height_tray*N_real) % (eqn. 6.13) *** question: N_theory or N_real?

% CALCULATING DIAMETER
        A_An_ratio = 0.8; % A_n/A is the fraction of the total area available for flow
        flood_frac = 0.6; % fraction of flooding velocity desired in the design
        
        R = 8.3145*10^-5;
        rho_l = 901% mass density of liquid
        rho_v = M_v * P / T / R % mass density of vapor
        
area_column = M_v/sqrt(rho_l*rho_v)/(flood_frac)*(A_An_ratio)*V; % [m^2]
diameter_column = 2*sqrt(area_column/pi) % [m]
% --------------------

%%
%comment

% CALCULATE HEAT EXCHANGER AREAS
% see (table 6.2) for U
% --------------------
% U is from (table 6.2) from Doherty's text
U = 350;

% changeT_A is the temperature difference between the two streams at end A,
% and changeT_B is the temperature difference between
%LMTD = ((TA-tB)-(TB-tA))/log((TA-tB)/(TB-tA));

Thot_in =100
Thot_out =90
Tcold_in = 30
Tcold_out = 90
% 
 LMTD = ((Thot_in - Tcold_out)-(Thot_out-Tcold_in))...
%     /log((Thot_in - Tcold_out)/(Thot_out-Tcold_in));
 Q = Q_r+Q_c; % **** is this correct?
LMTD = 70
A = Q/(U*LMTD)/3600; % **** check units of everything EVENTUALLY PUT IN METERS

% TO BE CONSISTANT WITH COSTING FUNCTION BELOW

 [ installed_cost_column, purchased_cost_heatex, installed_cost_heatex ] = EquipCosts_func( height_column, diameter_column, A );
% --------------------


% ********** need to
%       1) T hot and T cold
%       2) change U