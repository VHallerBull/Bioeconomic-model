close all
clear all

T=50*365;
% Setting the runs with respect to dispersal
alpha=7800;
c_dist=20000;
q=[.03/365,.007/365,.02/365];
% #############################
% % reserve_1 = 1 means no reserve, =2 means reserve in patch 1
% % Wage_Var = 1 means steady wage, =2 means variation based on equation,
% %      =3 means variation based on equation and randomness between individuals
% % Cost_Var = 1 means steady cost with distance, =2 randomly varied
% %      between individuals
% % Price_Var = 1 means steady price, =2 means variation based on equation,
% %      =3 means variation based on equation and randomness between individuals
% % q_Var = 1 means catchability fixed, =2 randomly varied between individuals

Wage_Var=1;
Cost_Var=1;
Price_Var=1;
q_Var=1;


%reserve=[5:50];

folder =['Trials_NoSassi2_workforceNonLinearFinallimit_Amanda_PriceFactor5'];
mkdir(folder)

reserve_names=[{'Large1'},{'Large2'},{'Large3'},{'Medium1'},{'Medium2'},'Medium3',...
    {'Small1'},{'Small2'},{'Small3'},{'TwoSmall1'},{'TwoSmall2'},{'TwoSmall3'},{'2Small3'}];
reserve_scenario=struct(reserve_names{1},[1:50],...
    reserve_names{2},[50:100],reserve_names{3},[100:150],...
    reserve_names{4},[10:40],reserve_names{5},[60:90],...
    reserve_names{6},[110:140],reserve_names{7},[20:30],...
    reserve_names{8},[70:80],reserve_names{9},[120:130],...
    reserve_names{10},[20:30,70:80],reserve_names{11},[20:30,120:130],...
    reserve_names{12},[70:80,120:130]);
Sassi=ones(1,T+1); %reserve is implemented at all times, change values o zero for when the reserve is out of effect
% for years=1:50
%     Sassi(years*350:years*350+14)=zeros(1,15);
% end
PriceFactor=5;

for i=2:2
%     if i==2
%         folder =['Trials_Sassib2'];
%         mkdir(folder)
%     end
    scenario=[2,1,1,1;...
        2,1,2,1];
    
    %         scenario=[1,1,1,1;...
    %             1,1,2,1;...
    %             1,1,3,1;...
    %             2,1,1,1;...
    %             3,1,1,1;...
    %             2,1,2,1;...
    %             3,1,3,1;...
    %             3,2,3,2];
    Wage_Var=scenario(i,1);
    Cost_Var=scenario(i,2);
    Price_Var=scenario(i,3);
    q_Var=scenario(i,4);
    
    
    entrants=3500;
    
    % x_out is T x 13:
    % Total dynamic welfare
    % harvest_sum by patch (3 vectors)
    % stock_sum (3)
    % no_go_sum (choice of outside alternative)
    % effort_sum trips by patch (3)
    % rent total (note difference from welfare)
    % avg pv rent (total pv rent divided equally over T for to conform
    
    % % reserve_1 = 1 means no reserve, =2 means reserve in patch 1
    
    tic
    reserve=0;
    [XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=PNASloop_nofigs_jns2(entrants, T, alpha,Wage_Var,Cost_Var, Price_Var,...
        q_Var,c_dist,q,reserve,Sassi,PriceFactor);
    Outcome(1,i)=struct('Biomass',XB,'Number',XN,'Harvest',H,'Effort',...
        effort,'Wage',alpha_av,'TotalIncome',tot_welfare,'Price',...
        price,'FishingMortality',F,'NumberNotFishing',...
        labour);
    TimeTaken(1)=toc;
%     
    for m=2:13
        m
       tic
        reserve=reserve_scenario.(reserve_names{m-1});
        [XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=PNASloop_nofigs_jns2(entrants, T, alpha,Wage_Var,Cost_Var, Price_Var,...
            q_Var,c_dist,q,reserve,Sassi,PriceFactor);
        Outcome(m+1,i)=struct('Biomass',XB,'Number',XN,'Harvest',H,'Effort',...
            effort,'Wage',alpha_av,'TotalIncome',tot_welfare,'Price',...
            price,'FishingMortality',F,'NumberNotFishing',...
            labour);
        TimeTaken(m)=toc;
        save([folder,'/Dataset3.mat'],'Outcome','-v7.3')
%         biomass_Total=[sum(Outcome(3,2).Biomass.Group1,1)',sum(Outcome(3,2).Biomass.Group2,1)',sum(Outcome(3,2).Biomass.Group3,1)'];
%         mortality_Total=[sum(Outcome(3,2).FishingMortality.Group1,2),sum(Outcome(3,2).FishingMortality.Group2,2),sum(Outcome(3,2).FishingMortality.Group3,2)];
%         for i=1:10950
%             Price_diff(i)=(Outcome(3,2).Price(i+1,1)-Outcome(3,2).Price(i,1))/Outcome(3,2).Price(i,1);
%             Harvest_diff(i)=(Outcome(3,2).Harvest(i+1,1)-Outcome(3,2).Harvest(i,1))/Outcome(3,2).Harvest(i,1);
%         end
        
    end
end
save([folder,'/Dataset3.mat'],'Outcome','TimeTaken','-v7.3')

%
% figure
% plot(1:size(XB.Group1,2),XB.Group1(1,:),1:size(XB.Group1,2),XB.Group2(1,:),1:size(XB.Group1,2),XB.Group3(1,:))
% legend('Group1','Group2','Group3')
% title('Biomass location1')
% saveas(gcf,[folder,'/Biomass.png'])
% 
% figure
% plot(1:size(XB.Group1,2),XN.Group1(1,:),1:size(XB.Group1,2),XN.Group2(1,:),1:size(XB.Group1,2),XN.Group3(1,:))
% legend('Group1','Group2','Group3')
% title('Number location1')
% saveas(gcf,[folder,'/Number.png'])
% 
% % figure
% % plot(1:size(H.Group1,2),H.Group1(1,:),1:size(H.Group1,2),H.Group2(1,:),1:size(H.Group1,2),H.Group3(1,:))
% % 
% % figure
% % plot(1:size(H.Group1,2),H.Group1(1,:)./XB.Group1(1,722:size(XB.Group1,2)),1:size(H.Group1,2),H.Group2(1,:)./XB.Group2(1,722:size(XB.Group1,2)),1:size(H.Group1,2),H.Group3(1,:)./XB.Group3(1,722:size(XB.Group1,2)))
% % 
% % figure
% % plot(1:size(H.Group1,2),H.Group1(2,:)./XB.Group1(2,722:size(XB.Group1,2)),1:size(H.Group1,2),H.Group2(2,:)./XB.Group2(2,722:size(XB.Group1,2)),1:size(H.Group1,2),H.Group3(2,:)./XB.Group3(2,722:size(XB.Group1,2)))
% 
% figure
% plot(1:size(price,1),price(:,1),1:size(price,1),price(:,2),1:size(price,1),price(:,3))
% legend('Group1','Group2','Group3')
% title('Price')
% saveas(gcf,[folder,'/Price.png'])
% 
% figure
% plot(1:size(price,1),H(:,1),1:size(price,1),H(:,2),1:size(price,1),H(:,3))
% legend('Group1','Group2','Group3')
% title('Harvest across all locations')
% saveas(gcf,[folder,'/Harvest.png'])
% 
% figure
% plot(1:size(effort,1),effort(:,1))
% title('labour')
% saveas(gcf,[folder,'/Labour.png'])
% 
% figure
% plot(1:size(effort,1),alpha_av(:,1))
% title('wage')
% saveas(gcf,[folder,'/Wage.png'])
% 
% figure
% Coordinates= csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/West_fromGIS.csv",1,2);
% Size=10*XN.Group1(:,end)./XN.Group1(:,1);
% Colour=10*XB.Group1(:,end)./XB.Group1(:,1);
% scatter(Coordinates(:,3),Coordinates(:,4),Size ,Colour,'fill')
% colorbar
% title('Group1')
% saveas(gcf,[folder,'/Spatial_Group1.png'])
% 
% figure
% Coordinates= csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/West_fromGIS.csv",1,2);
% Size=10*XN.Group2(:,end)./XN.Group2(:,1);
% Colour=10*XB.Group2(:,end)./XB.Group2(:,1);
% scatter(Coordinates(:,3),Coordinates(:,4),Size ,Colour,'fill')
% colorbar
% title('Group2')
% saveas(gcf,[folder,'/Spatial_Group2.png'])
% 
% figure
% Coordinates= csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/West_fromGIS.csv",1,2);
% Size=10*XN.Group3(:,end)./XN.Group3(:,1);
% Colour=10*XB.Group3(:,end)./XB.Group3(:,1);
% scatter(Coordinates(:,3),Coordinates(:,4),Size ,Colour,'fill')
% colorbar
% title('Group3')
% saveas(gcf,[folder,'/Spatial_Group3.png'])
% 
% figure
% scatter(labour(1:end-1,1),alpha_av(2:end,1))
% saveas(gcf,[folder,'/WageLabour.png'])
% 
% figure
% subplot(3,1,1)
% scatter(H(:,1)/H(1,1),price(:,1))
% subplot(3,1,2)
% scatter(H(:,2)/H(1,2),price(:,2))
% subplot(3,1,3)
% scatter(H(:,3)/H(1,3),price(:,3))
% saveas(gcf,[folder,'/HarvestPrice.png'])
% 
% figure
% plot(1:T+1,tot_welfare)
% saveas(gcf,[folder,'/Welfare.png'])
% 
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group1b_CalDataB.csv",XB.Group1(:,14601:15321));
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group2b_CalDataB.csv",XB.Group2(:,14601:15321));
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group3b_CalDataB.csv",XB.Group3(:,14601:15321));
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group1b_CalDataN.csv",XN.Group1(:,14601:15321));
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group2b_CalDataN.csv",XN.Group2(:,14601:15321));
% % csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group3b_CalDataN.csv",XN.Group3(:,14601:15321));
% 
% 
% 
% 
% 
% 
%save([folder,'/Basic.mat'],'X_out','X_out_r');

% fishing_mort_G1=mean((H.Group1./XB.Group1(:,722:4371))*365,2);
% fishing_mort_G2=mean(H.Group2./XB.Group2(:,722:4371)*365,2);
% fishing_mort_G3=mean(H.Group3./XB.Group3(:,722:4371)*365,2);
% csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Fishing_G1.csv",fishing_mort_G1);
% csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Fishing_G2.csv",fishing_mort_G2);
% csvwrite("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Fishing_G3.csv",fishing_mort_G3);




function [XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=PNASloop_nofigs_jns2(entrants, T, alpha,Wage_Var,Cost_Var, Price_Var,...
q_Var,c_dist,q,reserve,Sassi,PriceFactor)

discount=.07; % Discount rate
sigmult2=7000; % controls the variance of the Type I draws
logitshare=0; % 0 means use Type I draws, 1 means use logit probabilities (smoother)
n_locations=161;

% #######################
% BIOLOGICAL PARAMETERS (set in biological model code
% #######################

% #######################
% ECONOMIC PARAMETERS
% #######################
n=entrants; % Number of fishing vessels

% cathability
%q=.01; % can add more parameters later if interesting
% WE WILL TAKE q to be the mean of the catchability distribution
if q_Var ==2
    temp=rand(n,1);
else
    temp=.5*ones(n,1);% TURNS OFF HETEROGENEITY
end
q_ind=(q +(q).*(((2*temp).^0.5)-1)).*(temp<=0.5)+(q +(q).*(1-((2*(1-temp)).^0.5))).*(temp>0.5);

%figure()
%hist(q_ind)
%title 'Catchability Heterogeneity'


% prices and costs - assumes homogeneous harvest sector at beginning
%alpha=1; % outside opportunity
% WE WILL TAKE alpha to be the mean of the outside opp distribution
if Wage_Var==3
    temp=rand(n,1);
else
    temp=.5*ones(n,1);% TURNS OFF HETEROGENEITY
end
alpha_ind=(alpha +(alpha)*(((2*temp).^0.5)-1)).*(temp<=0.5)+(alpha +(alpha)*(1-((2*(1-temp)).^0.5))).*(temp>0.5);
%figure()
%hist(alpha_ind)
%title 'Opportunity Cost Heterogeneity'

base_p=[12900, 19200, 15300]; % base price
sigmult=1600; %standard deviation of price
%sigmult=0; % TURNS OFF PRICE VARIATION
c=0.2*45600;
if Cost_Var==2
    temp=rand(n,1);
else
    temp=.5*ones(n,1);% TURNS OFF HETEROGENEITY
end
c_dist_ind=(c_dist +(c_dist)*(((2*temp).^0.5)-1)).*(temp<=0.5)+(c_dist +(c_dist)*(1-((2*(1-temp)).^0.5))).*(temp>0.5);

%figure()
%hist(c_dist_ind)
%title 'Travel Cost Heterogeneity'


dist=repelem(1,n_locations); % Mean of the distance distribution for each site


% TURN ON DISTANCE HETEROGENEITY
dist_ind= (2.*rand(n,1)*dist)';

% TURN OFF DISTANCE HETEROGENEITY
%dist_ind= (repelem(1,n)'*dist)';

% #######################
% INITIAL CONDITIONS
% #######################


%Initial conditions (stock of fish in each patch  X(0)j)
% a percentage of carrying capacity

%Carrying capacity biomass (site*fish groups)
%Bmax=[repelem(14600/n_locations,n_locations)',repelem(42888/n_locations,n_locations)',repelem(365/n_locations,n_locations)'];
%fishing_mort= repelem(0,n_locations,3);
% Group1_B= repelem(30000/n_locations,n_locations)';
% Group1_N= repelem(840000/n_locations,n_locations)';
% Group2_B= repelem(85000/n_locations,n_locations)';
% Group2_N= repelem(8000000/n_locations,n_locations)';
% Group3_B= repelem(7300/n_locations,n_locations)';
% Group3_N= repelem(270000/n_locations,n_locations)';
Bmax=[35000/n_locations,150000*4/n_locations,10000*5/n_locations];

distance=csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/West_withDistance.csv",1,2);

% %initial stock biomass - run ode to find initial biomass
% for i=1:20*365
%     [t1,x1]= ode45(@biology1,[0 1/365],[Group1_B(:,end);Group1_N(:,end)],[],Group1_B,Bmax(:,1),Group1_N,distance,fishing_mort(:,1)');
%     %[t2,x2]= ode45(@biology2,[0 1/365],[Group2_B(:,end);Group2_N(:,end)],[],Group2_B,Bmax(:,1),Group2_N,distance,fishing_mort(:,2)');
%     %[t3,x3]= ode46(@biology3,[0 1/365],[Group3_B(:,end);Group3_N(:,end)],[],Group3_B,Bmax(:,1),Group3_N,distance,fishing_mort(:,3)');
%     Group1_B=[Group1_B,x1(end,1:n_locations)'];
%     Group1_N=[Group1_N,x1(end,n_locations+1:end)'];
%     %Group2_B=[Group2_B,x2(end,1:n_locations)'];
%     %Group2_N=[Group2_N,x2(end,n_locations+1:end)'];
%     %Group3_B=[Group3_B,x3(end,1:n_locations)'];
%     %Group3_N=[Group3_N,x3(end,n_locations+1:end)'];
% end

Group1_CalDataB = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group1b_CalDataB.csv");
Group1_CalDataN = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group1b_CalDataN.csv");
Group2_CalDataB = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group2b_CalDataB.csv");
Group2_CalDataN = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group2b_CalDataN.csv");
Group3_CalDataB = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group3b_CalDataB.csv");
Group3_CalDataN = 0.9*csvread("C:/Users/uqvhall1/OneDrive - The University of Queensland/documents/Bioeconomic model/data/Group3b_CalDataN.csv");

XB=struct('Group1',Group1_CalDataB,'Group2',Group2_CalDataB,'Group3',Group3_CalDataB);
XN=struct('Group1',Group1_CalDataN,'Group2',Group2_CalDataN,'Group3',Group3_CalDataN);
j=7;
 

% #############################
% % reserve_1 = 1 means no reserve, =2 means reserve in patch 1
%reserve_1=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function does all of the calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=pnasmarinereserve2(T,discount,...
    sigmult,n,n_locations,q_ind,alpha_ind,c_dist_ind, distance,dist,dist_ind,XB,XN,Bmax,base_p,c,...
    sigmult2,logitshare,Wage_Var,Cost_Var, Price_Var,q_Var,reserve,Sassi,PriceFactor);

%no_go_sum=n-sum(effort_sum')';

t=1:T;


%stock_tot=stock_sum(:,1)+ stock_sum(:,2)+ stock_sum(:,3);
%effort_tot=effort_sum(:,1)+ effort_sum(:,2)+ effort_sum(:,3);
%harvest_tot=harvest_sum(:,1)+harvest_sum(:,2)+harvest_sum(:,3);
%rent_tot=rent_sum(:,1)+rent_sum(:,2)+rent_sum(:,3);


%avgpvrentsum=sum(exp(-discount*t).*rent_tot')*ones(T,1);
%avg_dyn_welfare=tot_dyn_welfare/n;
%x_out=[tot_dyn_welfare', harvest_sum,stock_sum,no_go_sum,effort_sum,rent_tot,avgpvrentsum,price',labour',alpha_av'];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function file that runs the effort and population dynamics
% PNAS Marine reserve paper adjusted to simutaneously adjus price and
% harvest
% April 13, 2009
% Marty edits April 23, 2009
% Vanessa edits March 12, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [XB,XN,harvest,effort,alpha_av,tot_welfare,price,F,labour,Fit]=pnasmarinereserve2(T,discount,...
    sigmult,n,n_locations,q_ind,alpha_ind,c_dist_ind, distance,dist,dist_ind,XB,XN,Bmax,base_p,c,sigmult2,logitshare,Wage_Var,Cost_Var, Price_Var,...
q_Var,reserve,Sassi,PriceFactor);

% Calculating the price series
for i=2:T+1
    i;
    
    
    
    %Set starting values for price and harvest
    if i==2
        PT=base_p; %+ sigmult*randn();
    else
        PT=PT1; %+ sigmult*randn();
    end
    
    if PT<=0
        PT=[0,0,0];
        warning("Price drops below zero")
    end
    
    
    if i==2
        WT=alpha_ind(1);
    else
        WT=WT1;
    end
    if WT<=0
        WT=1;
        warning("Wage drops below zero")
    end
    
    [HT,LT,XB,XN,H,E_all,welfare,fishing_mort]=discreteChoicePriceWage(WT,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi );
    
    if LT<=0
        LT=0;
        warning("Labour drops below zero")
    end
    
    if HT<=0
        HT=[0,0,0];
        warning("Harvest drops below zero")
    end
    
    %'Total Welfare Cost Reserve in Patch 1'
    tot_welfare(i)=welfare;
    
    if i==2
     alpha_av(1:2,1) =[WT;WT];
     price(1:2,1:3)=[PT;PT];
     labour(1:2,1)=[LT;LT];
     harvest(1:2,1:3)=[HT;HT];
    else
        alpha_av(i,1) =[WT];
        price(i,1:3)=[PT];
        labour(i,1)=[LT];
        harvest(i,1:3)=[HT];
    end
    fishing1(i,1:n_locations)=fishing_mort(1,:);
    fishing2(i,1:n_locations)=fishing_mort(2,:);
    fishing3(i,1:n_locations)=fishing_mort(3,:);
    
    %discrete choice used to be here
    
    %tot_dyn_welfare(i)=sum(dyn_welfare);
    
    
    %effort_sum(i,1:n_locations)=E;

    
    harvest_G1_sum(i,1:n_locations)=sum(H.Group1,1);
    harvest_G2_sum(i,1:n_locations)=sum(H.Group2,1);
    harvest_G3_sum(i,1:n_locations)=sum(H.Group3,1);
    
    %rent_sum(i,1:n_locations)=sum(PT*H-(c+c_dist_ind.*dist1_ind +alpha_ind).*E);
    effort(i,1:n_locations+1)=sum(E_all,1);
    
    if Wage_Var==1
        if Price_Var==1
            PT1=PT;
            HT1=HT;
            LT1=LT;
            WT1=WT;
        else
            harvest1=price(i,1)-1.33*price(i,1)*(harvest(i,1)-harvest(i-1,1))/harvest(i-1,1);
            harvest2=price(i,2)-1.33*price(i,2)*(harvest(i,2)-harvest(i-1,2))/harvest(i-1,2);
            harvest3=price(i,3)-1.33*price(i,3)*(harvest(i,3)-harvest(i-1,3))/harvest(i-1,3);
            
            PT1=[harvest1,harvest2,harvest3];
            HT1=HT;
            LT1=LT;
            WT1=WT;
        end
    else if Price_Var==1
            WT1=alpha_av(i,1)-1.13*alpha_av(i,1)*(labour(i,1)-labour(i-1,1))/(labour(i-1,1)+50000);
            PT1=PT;
            LT1=LT;
            HT1=HT;
        else
             Harvest_diff1=(harvest(i,1)-harvest(i-1,1))/((harvest(i-1,1)));%+harvest(i-1,1))/2);
             Harvest_diff2=(harvest(i,2)-harvest(i-1,2))/((harvest(i-1,2)));%+harvest(i-1,2))/2);
             Harvest_diff3=(harvest(i,3)-harvest(i-1,3))/((harvest(i-1,3)));%+harvest(i-1,3))/2);
%              harvest1=price(i,1)*exp(-1.276*(Harvest_diff1));
%              harvest2=price(i,2)*exp(-1.276*(Harvest_diff2));
%              harvest3=price(i,3)*exp(-1.276*(Harvest_diff3)); 
             
%              harvest1=3.5965*price(i,1)*exp(-1.276*(Harvest_diff1+1));
%              harvest2=3.5965*price(i,2)*exp(-1.276*(Harvest_diff2+1));
%              harvest3=3.5965*price(i,3)*exp(-1.276*(Harvest_diff3+1));             
             if Harvest_diff1>1.4/3.2
                 Harvest_diff1=1.4/3.2;
             end
             if Harvest_diff2>1.4/3.2
                 Harvest_diff2=1.4/3.2;   
             end
              if Harvest_diff3>1.4/3.2
                  Harvest_diff3=1.4/3.2;
             end
             
             harvest1=price(i,1)-1.4*price(i,1)*Harvest_diff1+1.6*price(i,1)*Harvest_diff1^2;
             harvest2=price(i,2)-1.4*price(i,2)*Harvest_diff2+1.6*price(i,2)*Harvest_diff2^2;
             harvest3=price(i,3)-1.4*price(i,3)*Harvest_diff3+1.6*price(i,3)*Harvest_diff3^2;
             
             if harvest1>PriceFactor*12900
                 harvest1=PriceFactor*12900;
             end
             if harvest2>PriceFactor*19200
                 harvest2=PriceFactor*19200;
             end
             if harvest3>PriceFactor*15300
                 harvest3=PriceFactor*15300;
             end
 
             
             
%             harvest1=price(i,1)-1.33*price(i,1)*(harvest(i,1)-harvest(i-1,1))/(harvest(i-1,1));
%             harvest2=price(i,2)-1.33*price(i,2)*(harvest(i,2)-harvest(i-1,2))/(harvest(i-1,2));
%             harvest3=price(i,3)-1.33*price(i,3)*(harvest(i,3)-harvest(i-1,3))/(harvest(i-1,3));

             PT1=[harvest1,harvest2,harvest3];
             HT1=HT;
             WT1=alpha_av(i,1)-1.13*alpha_av(i,1)*(labour(i,1)-labour(i-1,1))/(labour(i-1,1)+50000);
             LT1=LT;
%            [x,fval]=fsolve(@(HP)FunPrice(HP,alpha_ind,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT,HT,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi),[HT,PT]);
%            PT1=[x(4),x(5),x(6)];
%            HT1=[x(1),x(2),x(3)];

        end
        
 
        
        
%         if X01<=0
%             X01=0;
%             warning("Location 1 extinct")
%         end
%         if X02<=0
%             X02=0;
%             warning("Location 2 extinct")
%         end
%         if X03<=0
%             X03=0;
%             warning("Location 3 extinct")
%         end
        
        
        
    end
    alpha_av(i,1:2) =[WT,WT1];
    price(i,1:6)=[PT,PT1];
    labour(i,1:2)=[LT,LT1];
    harvest(i,1:6)=[HT,HT1];
    Fit(i,1)=1;%:size(fval,2))=fval;
end

H2= struct('Group1',harvest_G1_sum','Group2',harvest_G2_sum','Group3',harvest_G3_sum');
F=struct('Group1',fishing1,'Group2',fishing2,'Group3',fishing3);

end

function F=FunPrice(HP,alpha_ind,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT,HT,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi)

% F(1)=PT(1)-HP(1)-1.4*PT(1)*(HP(4)-HT(1))/HT(1)+1.6*PT(1)*((HP(4)-HT(1))/HT(1))^2;
% F(2)=PT(2)-HP(2)-1.4*PT(2)*(HP(5)-HT(2))/HT(2)+1.6*PT(2)*((HP(5)-HT(2))/HT(2))^2;
% F(3)=PT(3)-HP(3)-1.4*PT(3)*(HP(6)-HT(3))/HT(3)+1.6*PT(3)*((HP(6)-HT(3))/HT(3))^2;

Harvest_diff1=(HP(4)-HT(1))/HT(1);
Harvest_diff2=(HP(5)-HT(2))/HT(2);
Harvest_diff3=(HP(6)-HT(3))/HT(3);
F(1)=3.5965*PT(1)*exp(-1.276*(Harvest_diff1+1));
F(2)=3.5965*PT(2)*exp(-1.276*(Harvest_diff2+1));
F(3)=3.5965*PT(3)*exp(-1.276*(Harvest_diff3+1));    

price=[HP(1),HP(2),HP(3)];
harvest=discreteChoicePriceWage(alpha_ind,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,price,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi);

F(4)=-HP(4)+harvest(1);
F(5)=-HP(5)+harvest(2);
F(6)=-HP(6)+harvest(3);

end

function F=FunWage(HP,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT,LT,WT,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi);

[HTI,LTI]=discreteChoicePriceWage(HP(2),q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi);
F(1)=-HP(1)+LTI;
if LT==0
    F(2)=WT-HP(2);
else F(2)=WT-HP(2)-1.13*WT*(HP(1)-LT)/LT;
end
end

function F=FunPriceWage(HP,PT,HT,LT,WT,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi);

price=[HP(4),HP(5),HP(6)];
[harvest,LTI]=discreteChoicePriceWage(HP(4),q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,price,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi);
F(1)=-HP(1)+harvest(1);
F(2)=-HP(2)+harvest(2);
F(3)=-HP(3)+harvest(3);
if HT(1)==0
    F(4)=PT(1)-HP(4);
else F(4)=PT(1)-HP(4)-1.33*PT(1)*(HP(1)-HT(1))/HT(1);
end
if HT(2)==0
    F(5)=PT(2)-HP(5);
else F(5)=PT(2)-HP(5)-1.33*PT(2)*(HP(2)-HT(2))/HT(2);
end
if HT(1)==0
    F(6)=PT(3)-HP(6);
else F(6)=PT(3)-HP(6)-1.33*PT(3)*(HP(3)-HT(3))/HT(3);
end
F(7)=-HP(7)+LTI;
if LT==0
    F(2)=WT-HP(2);
else F(2)=WT-HP(2)-1.13*WT*(HP(1)-LT)/LT;
end
end
  

%discrete choice model

function [HT1,LT1,XB,XN,H,E_all,dyn_welfare,fishing_mort]=discreteChoicePriceWage(WT1,q_ind,XB,XN,Bmax,c,c_dist_ind,dist_ind,PT1,n,distance,n_locations,sigmult2,logitshare,Wage_Var,Price_Var,reserve,Sassi )
 % Running the continuous time ordinary differential equations
 n=size(q_ind,1);
 
 if Wage_Var==3
     temp=rand(n,1);
 else
     temp=.5*ones(n,1);% TURNS OFF HETEROGENEITY
 end
 alpha_ind=(WT1 +(WT1).*(((2*temp).^0.5)-1)).*(temp<=0.5)+(WT1 +(WT1).*(1-((2.*(1-temp)).^0.5))).*(temp>0.5);
 
 if Price_Var==3
     sigmult=0;
 else
     sigmult=1000;% TURNS OFF randomness
 end
 
 p=PT1 + sigmult*randn();
 %p_ind=(PT1 +(PT1).*(((2.*temp).^0.5)-1)).*(temp<=0.5)+(PT1 +(PT1).*(1-((2.*(1-temp)).^0.5))).*(temp>0.5);

 
 X=[XB.Group1(:,end)';XB.Group2(:,end)';XB.Group3(:,end)'];
 
 
 % wage depending on choice
    vno=alpha_ind*8;
    vgo=q_ind.*p*X*1000-c-c_dist_ind.*dist_ind'; %p=1*groups; q_ind=entrants*1; X=location*groups;  
%     if any(reserve)
%         vgo(:,reserve)=repelem(0,size(vgo,1),size(reserve,2));
%     end
 
%     vno=alpha_ind*8;
%     vgo1=p_ind.*q_ind.*X01-c-c_dist_ind.*dist1_ind;
%     vgo2=p_ind.*q_ind.*X02-c-c_dist_ind.*dist2_ind;
%     vgo3=p_ind.*q_ind.*X03-c-c_dist_ind.*dist3_ind;

    %'Individual Welfare Cost Reserve in Patch 1'
    %welfare=log(exp(vno)+exp(vgo1)+exp(vgo2)+exp(vgo3))-log(exp(vno)+exp(vgo2)+exp(vgo3));
    welfare=1;
    
    sigmult2=8000;
    randmat=sigmult2*(-log(-log(rand(n,n_locations+1)))); % Type I Extreme Value Errors
    
   
    U0=vno+randmat(:,1); 
    U=vgo+randmat(:,2:end);
    
    if any(reserve)
        if Sassi(size(XN.Group1,2)-720)==1
            U(:,reserve)=repelem(0,size(vgo,1),size(reserve,2));
        end
    end

    E=zeros(n,n_locations);
    E_all=zeros(n,n_locations+1);

if logitshare==0; % means we use Type I Draws
    UALL=[U0,U];
    [maxU,ychoice]=max(UALL');
    
    LinIdx = sub2ind(size(E_all),[1:n]', ychoice'); 
    E_all(LinIdx)=1; % zeros besides location where they fish
    E0no=E_all(:,1);
    E=E_all(:,2:n_locations+1);
    dyn_welfare= sum([vno,vgo].*E_all,'all'); %log(exp(vno)+ exp(vgo1)+ exp(vgo2)+exp(vgo3));
else  % Means we use shares from the logit probabilities (smoother)
        UALL=[U0,U];
        [maxU,ychoice]=max(UALL');
        vno_exp=exp(vno);
        vgo_exp=exp(vgo);
        E0no=(vno_exp)./(vno_exp+sum(vgo_exp,2));
        E= (vgo_exp)./(vno_exp+sum(vgo_exp,2));

        dyn_welfare=log(sum(vno_exp)+sum(vgo_exp,'all')); 
end
    
    Group1_H=q_ind(:,1).*E*X(1,:)';
    Group2_H=q_ind(:,2).*E*X(2,:)';
    Group3_H=q_ind(:,3).*E*X(3,:)';
    H=struct('Group1',Group1_H,'Group2',Group2_H,'Group3',Group3_H);
    HT1=[sum(Group1_H,'all'),sum(Group2_H,'all'),sum(Group3_H,'all')];
    % H_loc=[sum(Group1_H,1),sum(Group2_H,1),sum(Group3_H,1)];
    LT1=sum(E0no);
    
    fishing_mort=[sum(q_ind(:,1)*365.*E,1);sum(q_ind(:,2)*365.*E,1);sum(q_ind(:,3)*365.*E,1)];
    
    [t1,x1]= ode45(@biology1,[0 1/365],[XB.Group1(:,end);XN.Group1(:,end);],[],XB.Group1,repelem(Bmax(1,1),161)',XN.Group1,distance,fishing_mort(1,:));
    [t2,x2]= ode45(@biology2,[0 1/365],[XB.Group2(:,end);XN.Group2(:,end);],[],XB.Group2,repelem(Bmax(1,2),161)',XN.Group2,distance,fishing_mort(2,:));
    [t3,x3]= ode45(@biology3,[0 1/365],[XB.Group3(:,end);XN.Group3(:,end);],[],XB.Group3,repelem(Bmax(1,3),161)',XN.Group3,distance,fishing_mort(3,:));
    
    
    XB= struct('Group1',[XB.Group1,x1(end,1:n_locations)'],'Group2',[XB.Group2,x2(end,1:n_locations)'],'Group3',[XB.Group3,x3(end,1:n_locations)']);
    XN= struct('Group1',[XN.Group1,x1(end,n_locations+1:end)'],'Group2',[XN.Group2,x2(end,n_locations+1:end)'],'Group3',[XN.Group3,x3(end,n_locations+1:end)']);

end



function dydt=biology1(t,X,B,Bmax,N,distance,fishing_mort)



days = size(B,2); %how many days have been run previously
timepoint = rem(days+20,30); %computes the remainder of next timestep over 365 days,i.e.it is zero once a year
habitat_avail=repelem(1,size(B,1));
habitat_quality=repelem(1,size(B,1));



%% Set biological parameters for the group
natural_mort = 0.3; %natural mortality
s = exp(-natural_mort); %proportion of natural adult survival per year
p = 0.5; %Brody growth coefficent (equivalent to K)
h = 0.7; %steepness of the recruitment function (defined as the percentage of maximum recruitment at 20% of the unfished biomass)
wP = 0.00012; %mean weight of young fish prior to recruitment
wR = 0.0002; %mean weight of finish at recruitment (reaching maturity or enetering fishery
w_inf = 0.009; %weight at age infinity

age = 550; %age in days from spawning to recruitment
sigma_A= 1.5; %distance that the species can travel based on home range (single
            %number)
sigma_L= 36;
Bt=X(1:161);



%%Caculate dispersal matrix
dispmat1= exp(-distance.^2./(2*sigma_L^2));
dispmat= dispmat1./repmat(sum(dispmat1,2),[1,size(Bt)]);


%%Recruitement
if timepoint==0 %Time at which recruitment occurs

   % Calculate reference points
    B0s = Bmax .* habitat_avail' .* habitat_quality';
    %B0 = nansum(B0s);
    S0s = nansum((B0s' .* dispmat),1); % unfished larval settlement
    %R0s = B0s.*exp(-natural_mort)/wR;
    R0s = (1 - s' - (p * s') + (p * s'.^2)) .* B0s ./ (wR - p * s' * wP);  % normalized virgin recruitment
    SPRcrit = (1 - h) / (4*h); % critical replacement threshold
    corrFac = max(eig(dispmat)); % correct SPRcrit if dispersal matrix is real or simulated (i.e. not idealized)
    compratCorr = 1 / (SPRcrit ./ corrFac); % apply correction to compensation ratio if needed
    alpha = compratCorr .* R0s' ./ S0s; % calculate Beverton-Holt parameter 1 (after Walters et al 2007)
    beta = (compratCorr - 1) ./ S0s; % calculate Beverton-Holt parameter 2 (after Walters et al 2007)
    R0s = alpha ./ (1 + beta .* S0s) .* S0s; % calculate virgin recruitment
    
    if age>=days
        Bpast= B(:,1);
    else
        Bpast= B(:,days-age); %set to ? days in the past
    end
   
    
    
    %%Calculate recruitement
    Sin = nansum((Bpast .* dispmat),1); % larval output *Check size of matrix/vectors
    Rpast = alpha ./ (1 + beta .* Sin) .* Sin; %recruitment
    Rpast=Rpast;
 
else Rpast=repelem(0,161);
end


Recruitment = wR.*Rpast; 

%%Set mortality
Zt = fishing_mort + repelem(natural_mort,size(fishing_mort,2));
Mortality = Zt'.*Bt;

%%Calculate adult movement
prop_moved= exp(-distance.^2./(2*sigma_A^2));%.*(habitat_avail'./Bt);
prop_moved_norm= prop_moved./sum(prop_moved,2);
Adults_moved = Bt.*prop_moved_norm;
if sigma_A==0
    Adults_moved=zeros(size(Adults_moved));
end

Emigration = sum(Adults_moved,2)-Adults_moved(:,1);
Immigration = sum(Adults_moved,1)'-Adults_moved(:,1);

%%Calculate growth of existing population
Nt=X(162:322)';
Growth = p.*w_inf.*Nt;     

%%Combine all and calculate biomass
 dydt = [Recruitment' + Growth' + Immigration - Emigration - Mortality;Rpast' - Zt'.*Nt'];
 %this creates a matrix with each row being the biomass of each cell,
 %followed by the number in each cell; the columns are time



end


function dydt=biology2(t,X,B,Bmax,N,distance,fishing_mort)


days = size(B,2); %how many days have been run previously
timepoint = rem(days+25,30); %computes the remainder of next timestep over 365 days,i.e.it is zero once a year

habitat_avail=repelem(1,size(B,1));
habitat_quality=repelem(1,size(B,1));


%% Set biological parameters for the group
natural_mort = 0.1; %natural mortality
s = exp(-natural_mort); %natural_mort; %proportion of natural adult survival per year
p = 0.25; %Brody growth coefficent (equivalent to K)
h = 0.5; %steepness of the recruitment function (defined as the percentage of maximum recruitment at 20% of the unfished biomass)
wP = 0.00012; %mean weight of young fish prior to recruitment
wR = 0.0002; %mean weight of fish at recruitment (reaching maturity or enetering fishery
w_inf = 0.005; %weight at age infinity

age = 470; %age in days from spawning to recruitment
sigma_A= 5; %distance that the species can travel based on home range (single
            %number)
sigma_L= 14;
Bt=X(1:161);



%%Caculate dispersal matrix
dispmat1= exp(-distance.^2./(2*sigma_L^2));
dispmat= dispmat1./repmat(sum(dispmat1,2),[1,size(Bt)]);


%%Recruitement
if timepoint==0 %Time at which recruitment occurs
    % Calculate reference points
    B0s = Bmax .* habitat_avail' .* habitat_quality';
    %B0 = nansum(B0s);
    S0s = nansum((B0s' .* dispmat),1); % unfished larval settlement
    %R0s = B0s.*natural_mort;
    R0s = (1 - s' - (p * s') + (p * s'.^2)) .* B0s ./ (wR - p * s' * wP);  % normalized virgin recruitment
    SPRcrit = (1 - h) / (4*h); % critical replacement threshold
    corrFac = max(eig(dispmat)); % correct SPRcrit if dispersal matrix is real or simulated (i.e. not idealized)
     compratCorr = 1 / (SPRcrit ./ corrFac); % apply correction to compensation ratio if needed
    alpha = compratCorr .* R0s' ./ S0s; % calculate Beverton-Holt parameter 1 (after Walters et al 2007)
    beta = (compratCorr - 1) ./ S0s; % calculate Beverton-Holt parameter 2 (after Walters et al 2007)
    R0s = alpha ./ (1 + beta .* S0s) .* S0s; % calculate virgin recruitment
    
    if age>=days
        Bpast= B(:,1);
    else
        Bpast= B(:,days-age); %set to ? days in the past
    end
   
    
    
    %%Calculate recruitement
    Sin = nansum((Bpast .* dispmat),1); % larval output *Check size of matrix/vectors
    Rpast = alpha ./ (1 + beta .* Sin) .* Sin; %recruitment
    Rpast=Rpast;
 
else Rpast=repelem(0,161);
end


Recruitment = wR.*Rpast; 



%%Set mortality
Zt = fishing_mort + repelem(natural_mort,size(fishing_mort,2));
Mortality = Zt'.*Bt;

%%Calculate adult movement
prop_moved= exp(-distance.^2./(2*sigma_A^2));%.*(habitat_avail'./Bt);
prop_moved_norm= prop_moved./sum(prop_moved,2);
Adults_moved = Bt.*prop_moved_norm;
if sigma_A==0
    Adults_moved=zeros(size(Adults_moved));
end

Emigration = sum(Adults_moved,2)-Adults_moved(:,1);
Immigration = sum(Adults_moved,1)'-Adults_moved(:,1);

%%Calculate growth of existing population
Nt=X(162:322)';
Growth = p.*w_inf.*Nt;     

%%Combine all and calculate biomass
 dydt = [Recruitment' + Growth' + Immigration - Emigration - Mortality;Rpast' - Zt'.*Nt'];
 %this creates a matrix with each row being the biomass of each cell,
 %followed by the number in each cell; the columns are time


end


function dydt=biology3(t,X,B,Bmax,N,distance,fishing_mort);



days = size(B,2); %how many days have been run previously
timepoint = rem(days+10,30); %computes the remainder of next timestep over 365 days,i.e.it is zero once a year

habitat_avail=repelem(1,size(B,1));
habitat_quality=repelem(1,size(B,1));



%% Set biological parameters for the group
natural_mort = 0.2; %natural mortality
s = exp(-natural_mort); %proportion of natural adult survival per year
p = 0.7; %Brody growth coefficent (equivalent to K)
h = 0.6; %steepness of the recruitment function (defined as the percentage of maximum recruitment at 20% of the unfished biomass)
wP = 0.00007; %mean weight of young fish prior to recruitment
wR = 0.0001; %mean weight of fish at recruitment (reaching maturity or enetering fishery
w_inf = 0.005; %weight at age infinity

age = 730; %age in days from spawning to recruitment
sigma_A= 50; %distance that the species can travel based on home range (single
            %number)
sigma_L= 15;
Bt=X(1:161);



%%Caculate dispersal matrix
dispmat1= exp(-distance.^2./(2*sigma_L^2));
dispmat= dispmat1./repmat(sum(dispmat1,2),[1,size(Bt)]);


%%Recruitement
if timepoint==0 %Time at which recruitment occurs
    
    % Calculate reference points
    B0s = Bmax .* habitat_avail' .* habitat_quality';
    %B0 = nansum(B0s);
    S0s = nansum((B0s' .* dispmat),1); % unfished larval settlement
    %R0s = B0s.*exp(-natural_mort)/wR;
    R0s = (1 - s' - (p * s') + (p * s'.^2)) .* B0s ./ (wR - p * s' * wP);  % normalized virgin recruitment
    SPRcrit = (1 - h) / (4*h); % critical replacement threshold
    corrFac = max(eig(dispmat)); % correct SPRcrit if dispersal matrix is real or simulated (i.e. not idealized)
    compratCorr = 1 / (SPRcrit ./ corrFac); % apply correction to compensation ratio if needed
    alpha = compratCorr .* R0s' ./ S0s; % calculate Beverton-Holt parameter 1 (after Walters et al 2007)
    beta = (compratCorr - 1) ./ S0s; % calculate Beverton-Holt parameter 2 (after Walters et al 2007)
    R0s = alpha ./ (1 + beta .* S0s) .* S0s; % calculate virgin recruitment
    
    if age>=days
        Bpast= B(:,1);
    else
        Bpast= B(:,days-age); %set to ? days in the past
    end
   
    
    
    %%Calculate recruitement
    Sin = nansum((Bpast .* dispmat),1); % larv bnvgbal output *Check size of matrix/vectors
    Rpast = alpha ./ (1 + beta .* Sin) .* Sin; %recruitment
    Rpast=Rpast;
 
else Rpast=repelem(0,161);
end


Recruitment = wR.*Rpast; 


%%Set mortality
Zt = fishing_mort + repelem(natural_mort,size(fishing_mort,2));
Mortality = (Zt'+p).*Bt;

%%Calculate adult movement
prop_moved= exp(-distance.^2./(2*sigma_A^2));%.*(habitat_avail'./Bt);
prop_moved_norm= prop_moved./sum(prop_moved,2);
Adults_moved = Bt.*prop_moved_norm;
if sigma_A==0
    Adults_moved=zeros(size(Adults_moved));
end

Emigration = sum(Adults_moved,2)-Adults_moved(:,1);
Immigration = sum(Adults_moved,1)'-Adults_moved(:,1);

%%Calculate growth of existing population
Nt=X(162:322)';
Growth = p.*w_inf.*Nt;     


%%Combine all and calculate biomass
 dydt = [Recruitment' + Growth' + Immigration - Emigration - Mortality;Rpast' - Zt'.*Nt'];
 %this creates a matrix with each row being the biomass of each cell,
 %followed by the number in each cell; the columns are time



end