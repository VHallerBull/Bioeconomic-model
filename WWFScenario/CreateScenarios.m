%Setup scenarios
i=1
for Reserve=1:9
    for DaysOpen=30:60
        for WeeksClosed=40:60
            ScenarioData(i,1)=Reserve;
            ScenarioData(i,2)=DaysOpen;
            ScenarioData(i,3)=WeeksClosed;
            i=i+1
        end
    end
end

csvwrite('ScenarioData.csv',ScenarioData)