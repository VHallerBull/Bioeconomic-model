%Plot WWF on HPC

%Overview for one reserve scenario across Sassi

%Download data from individual files
ReserveData= csvread("./Reserves.csv",1,0);
ScenarioData= csvread("./ScenarioData.csv");
ReserveNames= {'10_South','10_North','10_Central','20_South','20_North','20_Central','30_South','30_North','30_Central'};
i=1;
RangeStart=4558;
RangeEnd=5208;

for n=RangeStart:RangeEnd
    Reserve=ScenarioData(n,1);
    reserve_scenario= ReserveData(ReserveData(:,1)~=0,1)';
    reserve_name=ReserveNames{Reserve};
    Days_closed=ScenarioData(n,3)*7;
    Days_open=ScenarioData(n,2);
    PriceFactor=4;
    
    filename=[reserve_name,'_Closed',num2str(Days_closed),'_Open',num2str(Days_open),'_PriceFactor',num2str(PriceFactor),'.mat']; 
    if isfile(filename)
    Data=load(filename);
	Outcome=Data.Outcome;
    
    %Biomass year 1
    B_G1_Y1=(sum(Outcome.Biomass.Group1(:,365),1)-sum(Outcome.Biomass.Group1(:,1),1))/sum(Outcome.Biomass.Group1(:,1),1);
    B_G2_Y1=(sum(Outcome.Biomass.Group2(:,365),1)-sum(Outcome.Biomass.Group2(:,1),1))/sum(Outcome.Biomass.Group2(:,1),1);
    B_G3_Y1=(sum(Outcome.Biomass.Group3(:,365),1)-sum(Outcome.Biomass.Group3(:,1),1))/sum(Outcome.Biomass.Group3(:,1),1);
    
    %Biomass year 5
    B_G1_Y5=(sum(Outcome.Biomass.Group1(:,5*365),1)-sum(Outcome.Biomass.Group1(:,1),1))/sum(Outcome.Biomass.Group1(:,1),1);
    B_G2_Y5=(sum(Outcome.Biomass.Group2(:,5*365),1)-sum(Outcome.Biomass.Group2(:,1),1))/sum(Outcome.Biomass.Group2(:,1),1);
    B_G3_Y5=(sum(Outcome.Biomass.Group3(:,5*365),1)-sum(Outcome.Biomass.Group3(:,1),1))/sum(Outcome.Biomass.Group3(:,1),1);
    
    %Biomass year 10
    B_G1_Y10=(sum(Outcome.Biomass.Group1(:,10*365),1)-sum(Outcome.Biomass.Group1(:,1),1))/sum(Outcome.Biomass.Group1(:,1),1);
    B_G2_Y10=(sum(Outcome.Biomass.Group2(:,10*365),1)-sum(Outcome.Biomass.Group2(:,1),1))/sum(Outcome.Biomass.Group2(:,1),1);
    B_G3_Y10=(sum(Outcome.Biomass.Group3(:,10*365),1)-sum(Outcome.Biomass.Group3(:,1),1))/sum(Outcome.Biomass.Group3(:,1),1);
    
    %Harvest year 1
    H_G1_Y1=(Outcome.Harvest(365,1)-Outcome.Harvest(1,1))/Outcome.Harvest(1,1);
    H_G2_Y1=(Outcome.Harvest(365,2)-Outcome.Harvest(1,2))/Outcome.Harvest(1,2);
    H_G3_Y1=(Outcome.Harvest(365,3)-Outcome.Harvest(1,3))/Outcome.Harvest(1,3);
    
    %Harvest year 5
    H_G1_Y5=(Outcome.Harvest(5*365,1)-Outcome.Harvest(1,1))/Outcome.Harvest(1,1);
    H_G2_Y5=(Outcome.Harvest(5*365,2)-Outcome.Harvest(1,2))/Outcome.Harvest(1,2);
    H_G3_Y5=(Outcome.Harvest(5*365,3)-Outcome.Harvest(1,3))/Outcome.Harvest(1,3);
    
    %Harvest year 10
    H_G1_Y10=(Outcome.Harvest(3650,1)-Outcome.Harvest(1,1))/Outcome.Harvest(1,1);
    H_G2_Y10=(Outcome.Harvest(3650,2)-Outcome.Harvest(1,2))/Outcome.Harvest(1,2);
    H_G3_Y10=(Outcome.Harvest(3650,3)-Outcome.Harvest(1,3))/Outcome.Harvest(1,3);
    
    %Income 
    I_Y1=(Outcome.TotalIncome(365)-Outcome.TotalIncome(1))/Outcome.TotalIncome(1);
    I_Y5=(Outcome.TotalIncome(5*365)-Outcome.TotalIncome(1))/Outcome.TotalIncome(1);
    I_Y10=(Outcome.TotalIncome(3650)-Outcome.TotalIncome(1))/Outcome.TotalIncome(1);  
    
    Data(i,:)=[B_G1_Y1,B_G1_Y5,B_G1_Y10,B_G2_Y1,B_G2_Y5,B_G2_Y10,B_G3_Y1,...
        B_G3_Y5,B_G3_Y10,H_G1_Y1,H_G1_Y5,H_G1_Y10,H_G2_Y1,H_G2_Y5,...
        H_G2_Y10,H_G3_Y1,H_G3_Y5,H_G3_Y10,I_Y1,I_Y5,I_Y10];
    else Data(i,:)=zeros(1,21);
    end
    i=i+1
end

save([reserve_name,'_PlotData.mat'],'Data')

%Fig.1 1/5/10 year biomass change
figure
for m=1:9
subplot(3,3,m)
scatter(ScenarioData(RangeStart:RangeEnd,3),ScenarioData(RangeStart:RangeEnd,2),20,Data(:,m))
xlim([40 60])
ylim([30 60])
c=colorbar;
c.Label.String= 'Biomass(t)';
if m==1 | m==2 | m==3
    Group=1;
else if m==4 | m==5 | m==6
        Group=2;
    else Group=3;
    end
end
if m==1 | m==4 | m==7
    Year=1;
else if m==2 | m==5 | m==8
        Year =5;
    else Year=10;
    end
end

title(['Group ',num2str(Group),' at year ',num2str(Year)])
xlabel('Weeks closed')
ylabel('Days open')
end
saveas(gcf,[reserve_name,'_Figure1.png'])



%Fig.2 1/5/10 year harvest change
figure
for m=1:9
subplot(3,3,m)
scatter(ScenarioData(RangeStart:RangeEnd,3),ScenarioData(RangeStart:RangeEnd,2),40,Data(:,m+9))
xlim([40 60])
ylim([30 60])
c=colorbar;
c.Label.String= 'Harvest(t) per day';

if m==1 | m==2 | m==3
    Group=1;
else if m==4 | m==5 | m==6
        Group=2;
    else Group=3;
    end
end
if m==1 | m==4 | m==7
    Year=1;
else if m==2 | m==5 | m==8
        Year =5;
    else Year=10;
    end
end

title(['Group ',Group,' at year ',Year])
xlabel('Weeks closed')
ylabel('Days opened)
end
saveas(gcf,[reserve_name,'_Figure2.png'])

%Fig. 1/5/10 year income 
figure
for m=1:3
subplot(1,3,m)
scatter(ScenarioData(RangeStart:RangeEnd,3),ScenarioData(RangeStart:RangeEnd,2),40,Data(:,m+18))
xlim([40 60])
ylim([30 60])
c=colorbar;
c.Label.String= 'Income per day';

title(['Year ',m])
xlabel('Weeks closed')
ylabel('Days open')
end
saveas(gcf,[reserve_name,'_Figure3.png'])
