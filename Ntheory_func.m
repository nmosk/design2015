function [ N_theory_array, r_array ] = Untitled2( r_min, N_min )

r_max = (.95+r_min)/(1-.95); r_array =r_min:(r_max-r_min)/1000:r_max; % so that we go to nearly 0 to 1 on the x axis
N_theory_array = [];
y_array = [];

for n = 1:length(r_array)
syms N_theory_temp
N_theory_temp = solve( (N_theory_temp-N_min)/(N_theory_temp+1) == 0.75 .* (1 - ((r_array(n)-r_min)./(r_array(n)+1))^0.5688) , N_theory_temp);
y_value = (N_theory_temp-N_min)/(N_theory_temp+1);
N_theory_array(n) = N_theory_temp;
y_array(n) = y_value;
end
x_array = (r_array-r_min) ./ (r_array+1);

figure
hold on;
plot(x_array,y_array)
axis( [ 0, 1, 0, 1]);