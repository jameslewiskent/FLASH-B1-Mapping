%% Plots: First run Recon_Data.m
Localiser = dicomread('Three_Plane_Localiser.dcm');
Concat_Localiser = [flip(Localiser(:,:,1,1),2);Localiser(:,:,1,2);Localiser(:,:,1,3)];
figure('color','w');
imagesc(Concat_Localiser,[0 0.6*max(Localiser,[],'all')])
colormap('gray'); axis image off

load('SatTFL_Reconstructed_B1Maps_Images.mat');
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
load('Sandwich_Reconstructed_B1Maps_Images.mat');

figure('color','w'); imagesc(imtile(squeeze(FAperV_Combined(:,16,:) - FAperV(:,16,:,:))),[-0.2,0.2]);
axis image off
colorbar
colormap(bluewhitered)

figure('color','w'); tiledlayout(3,3,'TileSpacing','compact','Padding','none')
for Voltage_n = 1:size(Reference_Images,4)
nexttile; plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]),'.');
axis square
end

cmap = turbo(9);
figure('color','w'); tiledlayout(1,2); nexttile;
for Voltage_n = 1:size(Reference_Images,4)
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
lgd = legend(cellstr(num2str(Voltages')),'Location','Northwest','NumColumns',3); 
lgd.Title.String = 'Transmit Voltage [V]';
lgd.ItemTokenSize(1) = 10;
plot([0,150],[0,150],'k','HandleVisibility','off')
xlim([0 150]); ylim([0 150])
axis square
xlabel(['Reference Flip Angle, [',char(176),']'])
ylabel(['Measured Flip Angle, [',char(176),']'])

nexttile;
for Voltage_n = 1:size(Reference_Images,4)
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[])-reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
lgd = legend(cellstr(num2str(Voltages')),'Location','NorthEast','NumColumns',3); 
lgd.Title.String = 'Transmit Voltage [V]';
lgd.ItemTokenSize(1) = 10;
plot([0,150],[0,0],'k','HandleVisibility','off')
xlim([0 150]); ylim([-20 20])
axis square
xlabel(['Reference Flip Angle, [',char(176),']'])
ylabel(['Measured - Reference Flip Angle, [',char(176),']'])


%% Plot Lookup Tables
figure('color','w','Name','Lookup Tables'); 
load('SA2RAGE/sa2rage_lookup_table.mat');
plot(x_query,real(cos(x_query)),'k'); hold on
plot(x_query,fx_interp,'r');
load('Sandwich/sandwich_lookup_table.mat');
plot(x_query,fx_interp,'b'); hold on;
load('SatTFL/sattfl_lookup_table.mat');
plot(x_query,fx_interp,'k--'); hold on;

ylabel('Image Ratio, [a.u.]'); xlabel('Dynamic Range, [a.u.]')
ylim([-1 1])
legend('\it{Cosine}','SA2RAGE','Sandwich')
