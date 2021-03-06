function [ costs, costs_check ] = CostModel_func(N_real, V)
% Costs: column + heat exchangers + utilities
% Costs: Co*N^0.8*V^0.5 + C1*V^0.65+C2*V

C = [1*10^3 1*10^4 1*10^4];
costs = (C(1)*N_real.^0.8*V^0.5) + (C(2).*V.^0.65) + (C(3).*V); 
costs1 = (C(1)*N_real.^0.8*V.^0.5); % column cost
costs2 = (C(2)*V.^0.65); % heat exchanger cost
costs3 = (C(3)*V); % utilities cost
costs_check = [costs1 costs2 costs3];
end

