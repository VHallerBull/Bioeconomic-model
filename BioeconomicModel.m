function  [ScenarioName] = BioeconomicModel(T,reserve_scenario,reserve_name,Days_closed,Days_open,PriceFactor)

% #############################
% % Inputs
% % T=time to run model in days
% % reserve =  vector of all patches that have a reserve status   
% % Days_closed = number of days that the reserve is closed
% % Days_open = number of days that the reserve is opened to fishing
% %             following the closure
% % PriceFactor = multiple of teh original price that is set as the maximum
% %               price


% Setting the economic inputs
alpha=7800; %wage
c_dist=20000; %cost per distance
q=[.03/365,.007/365,.02/365]; %catchability coefficient for each species

% % Wage_Var = 1 means steady wage, =2 means variation based on equation,
% %      =3 means variation based on equation and randomness between individuals
% % Cost_Var = 1 means steady cost with distance, =2 randomly varied
% %      between individuals
% % Price_Var = 1 means steady price, =2 means variation based on equation,
% %      =3 means variation based on equation and randomness between individuals
% % q_Var = 1 means catchability fixed, =2 randomly varied between individuals

Wage_Var=2;
Cost_Var=1;
Price_Var=2;
q_Var=1;


folder =['Results/',reserve_name,'_Closed',num2str(Days_closed),'_Open',num2str(Days_open),'_PriceFactor',num2str(PriceFactor)];%folder to save results created through the function inputs
mkdir(folder)

%Set up Sassi reserve, 1=reserve open, 0=reserve closed
Start=1;
End=Days_closed+Days_open;
Sassi=[];
while length(Sassi)<T
    Sassi(Start:End)=zeros(1,Days_closed+Days_open);
    Sassi(Start+Days_closed:End)=ones(1,Days_open);
    Start=End+1;
    End=Start+Days_closed+Days_open-1;
end
Sassi=Sassi(1:T);
  
    
entrants=3500; %number of fishermen
    

reserve=reserve_scenario.(reserve_names{m-1});
[XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=PNASloop_nofigs_jns2(entrants, T, alpha,Wage_Var,Cost_Var, Price_Var,...
    q_Var,c_dist,q,reserve,Sassi,PriceFactor);
Outcome=struct('Biomass',XB,'Number',XN,'Harvest',H,'Effort',...
    effort,'Wage',alpha_av,'TotalIncome',tot_welfare,'Price',...
    price,'FishingMortality',F,'NumberNotFishing',...
    labour);
TimeTaken(m)=toc;
save([folder,'/Dataset_reserve.mat'],'Outcome','-v7.3')



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
Bmax=[35000/n_locations,150000*4/n_locations,10000*5/n_locations];

distance=csvread("./DataSelayar/West_withDistance.csv",1,2);

Group1_CalDataB = 0.9*csvread("./DataSelayar/Group1b_CalDataB.csv");
Group1_CalDataN = 0.9*csvread("./DataSelayar/Group1b_CalDataN.csv");
Group2_CalDataB = 0.9*csvread("./DataSelayar/Group2b_CalDataB.csv");
Group2_CalDataN = 0.9*csvread("./DataSelayar/Group2b_CalDataN.csv");
Group3_CalDataB = 0.9*csvread("./DataSelayar/Group3b_CalDataB.csv");
Group3_CalDataN = 0.9*csvread("./DataSelayar/Group3b_CalDataN.csv");

XB=struct('Group1',Group1_CalDataB,'Group2',Group2_CalDataB,'Group3',Group3_CalDataB);
XN=struct('Group1',Group1_CalDataN,'Group2',Group2_CalDataN,'Group3',Group3_CalDataN);
j=7;
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function does all of the calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XB,XN,H,effort,alpha_av,tot_welfare,price,F,labour,Fit]=pnasmarinereserve2(T,discount,...
    sigmult,n,n_locations,q_ind,alpha_ind,c_dist_ind, distance,dist,dist_ind,XB,XN,Bmax,base_p,c,...
    sigmult2,logitshare,Wage_Var,Cost_Var, Price_Var,q_Var,reserve,Sassi,PriceFactor);

%no_go_sum=n-sum(effort_sum')';

t=1:T;
         
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
    
    
    harvest_G1_sum(i,1:n_locations)=sum(H.Group1,1);
    harvest_G2_sum(i,1:n_locations)=sum(H.Group2,1);
    harvest_G3_sum(i,1:n_locations)=sum(H.Group3,1);
    

    effort(i,1:n_locations+1)=sum(E_all,1);
    
    if Wage_Var==1
        if Price_Var==1
            PT1=PT;
            HT1=HT;
            LT1=LT;
            WT1=WT;
        else
            Harvest_diff1=(harvest(i,1)-harvest(i-1,1))/((harvest(i-1,1)));
             Harvest_diff2=(harvest(i,2)-harvest(i-1,2))/((harvest(i-1,2)));
             Harvest_diff3=(harvest(i,3)-harvest(i-1,3))/((harvest(i-1,3)));
           
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
             Harvest_diff1=(harvest(i,1)-harvest(i-1,1))/((harvest(i-1,1)));
             Harvest_diff2=(harvest(i,2)-harvest(i-1,2))/((harvest(i-1,2)));
             Harvest_diff3=(harvest(i,3)-harvest(i-1,3))/((harvest(i-1,3)));
           
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
             
             PT1=[harvest1,harvest2,harvest3];
             HT1=HT;
             WT1=alpha_av(i,1)-1.13*alpha_av(i,1)*(labour(i,1)-labour(i-1,1))/(labour(i-1,1)+50000);
             LT1=LT;


        end       
        
        
    end
    alpha_av(i,1:2) =[WT,WT1];
    price(i,1:6)=[PT,PT1];
    labour(i,1:2)=[LT,LT1];
    harvest(i,1:6)=[HT,HT1];
    Fit(i,1)=1;
end

H2= struct('Group1',harvest_G1_sum','Group2',harvest_G2_sum','Group3',harvest_G3_sum');
F=struct('Group1',fishing1,'Group2',fishing2,'Group3',fishing3);

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
    Sin = nansum((Bpast .* dispmat),1); % larval output *Check size of matrix/vectors
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