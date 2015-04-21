% calculating N theoretical

% Function that calculates N_theory values for r_min to 1
% uses the Gilliland design method equation (eqn 4.56)
% also outputs a plot of (N-N_min)/(N+1) vs (r-r_min)/(r+1)

% function [ N_theory_array, r_array, N_theory ] = Ntheory_func( r_min, N_min, r )
function [ N_theory ] = Ntheory_func( r_min, N_min, r )
% r_max = (.95+r_min)/(1-.95); r_array =r_min:(r_max-r_min)/1000:r_max; % so that we go to nearly 0 to 1 on the x axis
% N_theory_array = [];
% y_array = [];
% 
% for n = 1:length(r_array)
% syms N_theory_temp
% N_theory_temp = solve( (N_theory_temp-N_min)/(N_theory_temp+1) == 0.75 .* (1 - ((r_array(n)-r_min)./(r_array(n)+1))^0.5688) , N_theory_temp);
% y_value = (N_theory_temp-N_min)/(N_theory_temp+1);
% N_theory_array(n) = N_theory_temp;
% y_array(n) = y_value;
% end
% x_array = (r_array-r_min) ./ (r_array+1);
syms N_theory
N_theory = solve( (N_theory-N_min)/(N_theory+1) == 0.75 .* (1 - ((r-r_min)./(r+1))^0.5688) , N_theory);

% figure
% hold on;
% plot(x_array,y_array)
% axis( [ 0, 1, 0, 1]);
end