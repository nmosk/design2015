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
F % feed molar flowrate [mol/hr]
q = 1; % fraction of feed that is liquid
z_F % composition across all phases
T % temperature in K
P % pressure in bar

% The species in distillate and bottoms
% 1 = species present ; 0 = species not present
species_D = [1 1 1 0 0 0 0 0]; 
species_B = [0 0 0 0 1 1 1 0];
% --------------------
V % vapor rate
zB
zD
%%

% CALCULATE D and B
% F = D + B
% basis of 100% recovery of LK in distillate and 100% recovery of JK in
% bottoms
% --------------------
D = F.*((zF(i)-zB(i))/(zD(i)-zB(i))); % (eqn 4.6) [mol/hr]
B = F.*((zD(i)-zF(i))/(zD(i)-zB(i))); % (eqn 4.7) [mol/hr]
% --------------------

% CALCULATE MINIMUM REFLUX AND REFLUX
% use Doherty's book, chapter 4 for pseudo 4-composition mixture if >4
% species
% --------------------

% calculating relative volatilities of multi-component mixture
% assuming that volatility is constant
% alpha = (y_LK/x_LK) / (y_HK/x_HK) = K_LK/K_HK

    function [RelVol] = RelVol_func(T); % where St is the reference component
    % for an AB/CD split, A = LK and C = HK
    % T is the LK and eB is the HK
    % where A=1 , B=2, C=3, D=4, etc
    % where A=Toluene , B=Benzene, C=Ethylbenzene, D=Styrene
    RelVol = [# # # 1]; % relative volatilities [A B C D] 
    % where D=1 since St is the reference component

r_min = ( (alpha(3)*x(1)/ (alpha(1)-alpha(3))) + ((alpha(3)*(x(2)+x(3)))/(alpha(2)-alpha(3))) ) / ((x(1)+x(2)) * (1+(x(1)*(x(3)+x(4)))))% minimum reflux for an AB/CD split

    % FOR BINARY -
    % % alpha = (y_i/x_i) / (y_j/x_j) = K_i/K_j
    % alpha = (y(1)/x(1)) / (y(2)/x(2)); 
    % % where i is the more volatile component and j is the less volatile
    % % component
    % r_min = (alpha(1)-1) * ( (x_D/x_F)-(alpha(1)*(1-x_D)/(1-x_F)) ); %
    % % Underwood equation for a binary mixture from "Distillation: Fundamentals
    % % and Principles" edited by Andrzej Gorak and Eva Sorensen (eqn 4.71)
    % ------------

r = r_min*1.5;
% --------------------

% *** after this point, no more choice in reboiler ratio

s = (D/B)*(r+q)-(1-q); % (eqn 3.35)

%%

% CALCULATE MINIMUM NUMBER OF STAGES
% use Fenske equation
% --------------------
N_min = log((zD(i)/zD(j))*(zB(j)/zB(i)))/log(alpha(i,j)); % (eqn 4.15)
% --------------------

% CALCULATE THEORETICAL AND REAL NUMBER OF STAGES
% use FUG method (p. 136)
% --------------------
N_theory
N_real = 2*N_theory;
% --------------------

% CALCULATE VAPOR RATES [mol/hr] IN COLUMN
% --------------------
v_B = s*B; % in bottoms
v_T = (r+1)*D; % in tops

% CROSS-CHECK: v_B-v_T = (q-1)*F (eqn 3.39)
%              when q = 1, v_B = v_T = V
% --------------------

% CALCULATE HEAD LOADS (saturated liquid products)
% lambda is the latent heat of vaporization
% estimate by taking an arithmatic average
% --------------------
[ lambda_D, lambda_B ] = HeatVap_func(x_D,x_B)

Q_c = lambda_D*v_T; % [J/hr]
Q_r = lambda_B*v_B; % [J/hr]
% --------------------

%%

% CALCULATE COLUMN DIMENSIONS
% height with (eqn. 6.13)
% diameter with (eqn. 6.12)
% --------------------


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

%   height_theory is approximately 30 to 60 cm (12 to 24 inches) but can be
%   smaller or larger for some columns
    height_theory = 0.5   % spacing between trays [m]

        M_v % molecular weight of the vapor
        M_l % molecular weight of the liquid
        c_o = 439; % [m/h] % from (table 6.1)
        height_column_min = 3*height_theory;
    
    height_column = height_column_min+(height_theory*N_real); % (eqn. 6.13) *** question: N_theory or N_real?

% CALCULATING DIAMETER
        A_An_ratio = 0.8; % A_n/A is the fraction of the total area available for flow
        flood_frac = 0.6; % fraction of flooding velocity desired in the design
        rho_l % mass density of liquid
        rho_v % mass density of vapor
        
area_column = M_v/sqrt(rho_l*rho_v)/(flood_frac)*(A_An_ratio)*V; % [m^2]
diameter_column = 2*sqrt(area_column/pi); % [m]
% --------------------

%%

% CALCULATE HEAT EXCHANGER AREAS
% see (table 6.2) for U
% --------------------
% changeT_A is the temperature difference between the two streams at end A,
% and changeT_B is the temperature difference between 
LMTD = 
A = Q/(U*LMTD);
% --------------------