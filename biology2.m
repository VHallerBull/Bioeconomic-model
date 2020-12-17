%Continuous time delay difference model
%by Vanessa Haller-Bull, updated December 2020


%function to solve the biology over a certain timestep
%[0 1/365] -> calculate the biomass from timstep 0 to 1/365 (one day)
%[XB.Group2(:,end);XN.Group2(:,end);] -> starting values for x
%[] -> options, kept clear in ths case
%XB.Group2,repelem(Bmax(1,2),161)',XN.Group2,distance,fishing_mort(2,:) ->
%function inputs


[t2,x2]= ode45(@biology2,[0 1/365],[XB.Group2(:,end);XN.Group2(:,end);],[],XB.Group2,repelem(Bmax(1,2),161)',XN.Group2,distance,fishing_mort(2,:));



function dydt=biology2(t,X,B,Bmax,N,distance,fishing_mort)

%Input:
%t = vector of time steps
%X = vector of biomass at each location at teh start
%B = matrix of biomass of fish per location (rows) nd timesteps in the past
%(columns)
%Bmax = vector of maximum biomass possible in each loaction
%N = matrix of number of fish per location (rows) nd timesteps in the past
%(columns)
%distance = matrix of distances between locations
%fishing_mort = vector of fishing mortality at each location



days = size(B,2); %how many days have been run previously
timepoint = rem(days+25,30); %computes the remainder of next timestep over 365 days,i.e.it is zero once a year

habitat_avail=repelem(1,size(B,1)); %habitat kept constant at this point
habitat_quality=repelem(1,size(B,1));%habitat kept constant at this point


%% Set biological parameters for the group 2 consisting of grouper, snapper and emperors
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
sigma_L= 14;%distance that the lrval can travel(single number)

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
