function [H_E,ROI_BT, reac, V_ft, D_fact ,WC_CF ,PO_CF ,  TCI, H, D, FC,TI, SU, WC, Profit_BT, Profit_AT, C_F, Cashflow_d, Bond_Fin, D_CF, NPV_0, NPV_proj,NPV_percent,Depreciation] = conceptual_econ_V6(V, WC, EPF,X)
close all
F_m=1; F_p=1; F_c=F_m*F_p;
FR=.04;
TR=.48;
ER=0.08; % enterprise rate
MAS=1600;
an3=0.00;
an2=0.00;
an1=0.50;
a0=0.50;
alphasv=0.03; %alpha salvage value
format long
Nconst=2; % Nconstructions
Nop=10; %Noperations
alphaSU=.1; 

IF=2.18; %installation factor for labor, foundations, supports, etc
CR=.06; %Construction rate -- chosen independently

%%
a=1/3; b=6; %[m]
D=transpose(a:abs((a-b)./(length(V)-1)):b); %diameter in [m]
% Changes diameter from meters to feet
D=D*3.28; % [ft]

% Changes V from meters^3 to feet^3
V_ft=35.3147.*V;

% Calculates corresponding height to chosen D and 
H=V_ft./(pi.*(D./2).^2); % [ft]
P_BT=EPF;%_fuel; %have form EP_fuel


%%

Bond_Fin=zeros(10,length(D));
ROI_BT=[]; TCI=[]; FC=[];
NPV_0=[]; NPV_proj=[]; NPV_percent=[]; Depreciation=[]; Profit_AT=[];C_F=[];D_CF=[];


         %Discount factors
         D_fact=zeros(10,1);
         D_fact(1)=1/(1+ER); %discount factors with enterprise rate at year 1
         for k=2:10
             D_fact(k)=D_fact(k-1)/(1+ER); % discount factors with enterprise rate at subsequent years
         end
         
%for i=1:length(V)
    %for j=1:length(D)
        
        %use 'i' if same value for every D in a certain V
        %use 'j' if different value for every D in a certain V
        
        %Purchase Cost of Base Equipment
       
        reac= MAS./280.*101.9.*(D.^1.066).*(H.^0.82 ); % purchasing cost of reactor
        A_c= 2*pi*(D./2).*H+2*pi*(D./2).^2;  % area of cylinder 
        H_E= + MAS/280*101.3*A_c.^0.65 ;% purchasing cost of heat exchanger
%        sep= (MAS./280)*4.7.*D.^1.55.*H*F_c;
%        sep=4E6;
sep=(2.*10.^6./X);
        PCBE=reac +H_E + sep;
        
        
        %Finding SU "Startup Cost"
        ISBL=PCBE.*(F_c+IF); %Installation cost
        
        FC=2.28*ISBL ;
        FC=FC';
        %FC=FC(1,:)
        SU=0.1.*FC; 
        %SU=SU(1,:)
        %WC: raw material cost for 2 months
        %         TI(i,j)=(2.5*PCBE*(F_c+IF))+WC(i);

        WC=FC.*.2;
        %WC=WC(1,:)
       
        TI=FC*(1+.2+.1);
        %TI=TI(1,:)
        %OUTPUTS
        FC_n3=FC*an3;
        FC_n2=FC*an2; 
        FC_n1=FC*an1; 
        FC_0=FC*a0;
        
        FC_n3=transpose(FC_n3);
        FC_n2=transpose(FC_n2) ;
        FC_n1=transpose(FC_n1);
        FC_0=transpose(FC_0);
      
        F_d=[(1+CR)^3 (1+CR)^2 (1+CR)^1 1 1 1]; %Discount factors for Y-3 to Y-0 then WC, SU
       
        Cashflow_d=[F_d(1).*FC_n3 F_d(2).*FC_n2 F_d(3).*FC_n1 FC_0.*F_d(4) WC'.*F_d(5) SU'.*F_d(6)];
        Cashflow_d=Cashflow_d';
        TCI=abs(sum(Cashflow_d));
       % TCI=TCI(1,:)
        
        ROI_BT=P_BT'./TI; %percent -- make sure to look at (i,j) cell
        %         TCI(i,j)=(an3*FC(j)*(1+CR))+(an2*FC(j)*(1+CR))+(an1*FC(j)*(1+CR))+(a0*FC(j))+WC(i)+SU(j); %total capital investment
        % b coefficients
        b_1=0.8;
        b_2=0.9;
        b_3=0.95;
        b_coeff=[b_1 b_2 b_3];
        % calcs profit before taxes for the first 3 years 
       
        %Profit_BT=zeros(3,length(V));
%         for k=1:length(V)
%             for z=1:length(b_coeff)
%             Profit_BT(z,k)=P_BT(k)*b_coeff(z);
%             end
%         end

        Profit_BT = transpose(b_coeff)*transpose(P_BT);
        
%         Profit_BT_temp=zeros(7,length(V));
%         for k=1:length(V)
%             for z=1:7
%             Profit_BT_temp(z,k)=P_BT(k);
%             end
%         end
%         Profit_BT=[Profit_BT; Profit_BT_temp]; % profit before taxes for 10 years; array is 10 by length(D)
        
        Profit_BT = [Profit_BT; repmat(transpose(P_BT),7,1)];



%          for k=1:10
%               Bond_Fin(k,j)=-FR*TCI(j);
%          end
%        Bond_Fin=Bond_Fin(1,:); %should be array of 10 by length(D) BUT all 10 values in each column are same so truncated

        Bond_Fin = -FR*TCI;

         %Depreciation
         Depreciation=-0.1.*FC.*(1+alphaSU); % depriciation allowed -0.1*FC*(1+alpha_Start_Up_Capital)
         %should be 10 by length(D) but all 10 rows are same value

         
         %Profit after taxes
%          for k=1:10
%          Profit_AT(k,j)=(Profit_BT(k,j)+Bond_Fin(1,j)+Depreciation(1,j))*(1-TR); %profit after tax
%          end
         Bond_Fin = repmat(Bond_Fin,10,1) ;
        %Depreciation=Depreciation'
         Depreciation = repmat(Depreciation,10,1) ; 
         
      
        
         Profit_AT = (Profit_BT + Bond_Fin + Depreciation)*(1-TR);
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        %Cash flow
%          for k=1:10
%          C_F(k,j)=Profit_AT(k,j)-Depreciation(j); % cash flow
%          end
       
        
        C_F = Profit_AT - Depreciation;
                   
%Discounted cash flow
         D_CF=(repmat(D_fact,1,10000)).*C_F; % discounted cash flows
         
     
       
         
         %NPV_0
         
            sum_D_CF=sum(D_CF);
         
       %%    
         SV=FC*alphasv; %alpha salvage value
         
         Pay_Off_TCI=-TCI; %pay off TCI
         
         Profit_AT_SV= SV*(1-TR); % profit after tax for salvage value
         
         SV_CF = D_fact(end)*Profit_AT_SV; % discounted cashflow for salvage value
         WC_CF= WC*D_fact(end);%  discounted cashflow for WC
         PO_CF= Pay_Off_TCI*(1.00); %discounted cashflow for WC
         
         NPV_0=sum_D_CF+SV_CF+WC_CF+ PO_CF;
        %%
         %NPV_proj
         NPV_proj=NPV_0/((1+ER)^Nconst);
         
         %NPV
         
         NPV_percent=(NPV_proj./TCI./(Nconst+Nop)).*100;
        

end



