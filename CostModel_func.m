function [ costs, costs_check ] = CostModel_func(N_real, V)
% Costs: column + heat exchangers + utilities
% Costs: Co*N^0.8*V^0.5 + C1*V^0.65+C2*V

C = [1*10^3 1*10^5 1*10^4];
costs = (C(1)*N_real^0.8*V^0.5) + (C(2)*V^0.65) + (C(3)*V); % where 10 is 10 years because the third term is $/yr
costs1 = (C(1)*N_real^0.8*V^0.5);
costs2 = (C(2)*V^0.65);
costs3 = (C(3)*V);
costs_check = [costs1 costs2 costs3];
end

