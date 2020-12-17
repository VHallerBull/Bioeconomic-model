% Create scenarios
Reserve_scenarios=zeros(48,100);

%30% reserve
Size=[48,24,16,12,8];
Number=[1,2,3,4,6];
Reserve_names=[];
scenario=1;
for m=1:5
    for i=1:20
        Reserve=[];
        Location=linspace(1,161,161);
        for n=1:Number(m)
            Consecutive=2;
            while min(Consecutive)~=1
                Random=randi([1,size(Location,2)-Size(m)],1);
                Reserve_try=linspace(Location(Random), Location(Random)+Size(m)-1, Size(m));
                Consecutive=diff(Reserve_try);
            end
            Reserve=[Reserve, Reserve_try];
            Location=setdiff(Location,Reserve);
        end
        Reserve_scenarios(1:48,scenario)=Reserve';
        Reserve_name=['Reserves',num2str(m),'_Random',num2str(i)];
        ReserveNames={Reserve_names,Reserve_name};
        scenario=scenario+1;
    end
end

save('Scenario_Independent1.mat','ReserveNames','Reserve_scenarios')

