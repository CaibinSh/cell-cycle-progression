%% determine endpoint cell cycle

% Right after live cell imaging, we perform EdU labeling and stain DNA with Hoechst. 
% The input Hoechst intensity and nuc_size are measured by Hoechst staining
% EdU means EdU instensity
% Output cycle: 1 (G1-phase); 2 (S-phase); 3 (G2-phase)
% see cycler1 function

cycle=cycler1(Hoechst,EdU,nuc_size);

%% plot cell cycle distribution by clusters

% plot cell cycle distribution for all cells in 1st column
subplot(4,3,4)
c = categorical(cycle,[1 2 3],{'G1','S','G2'})
histogram(c)
ylabel('cell counts')
box on 
title('all cells')
hold on,

% plot 1st level clusters in 2nd column and 2nd level clusters in 3rd
% column
for ii =1:2

    temp_cycle = cycle(subgroup_ID(:,1)==ii)
    
    for iii=0:2
        clear temp_cycle_subgroup
        if iii ==0
            temp_cycle_subgroup = temp_cycle;
        else
            temp_cycle_subgroup = cycle(subgroup_ID(:,2)==ii*10+iii);
        end
        
        subplot(4,3,ii*6-4+iii*iii)
        hold on, 
        c = categorical(temp_cycle_subgroup,[1 2 3],{'G1','S','G2'})
        histogram(c)
        title(['subgroup: ',num2str(ii*10+iii)])
        ylabel('cell counts')
        box on
    end
end

print(gcf,'-dpdf', '-noui',['plot/','10.CellCycles_by_clusters ','(',date,').pdf'])
