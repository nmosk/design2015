% Antoines eq for vapor pressure
% source: http://www.eng.auburn.edu/~drmills/mans486/Diffusion%20Tube/Antoine_coefficient_table.PDF

function vap_P = vap_P(temp,comp)

T = temp;
switch (comp)
    
    case 'st'
        A = 6.96;
        B = 1446;
        C = 209;
    case 'eb'
        A = 6.95;
        B = 1424;
        C = 213;
    case 'to'
        A = 6.95;
        B = 1344;
        C = 219;
    case 'be'
        A = 6.9;
        B = 1211;
        C = 220;
end

logP = A - B ./(T + C); % mmHg
vap_P = exp(logP) / 750; % bar




