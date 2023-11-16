%% (First run Recon_Data.m)
%% Plot Figure 1 Ground Truth Maps

Localiser = dicomread('Three_Plane_Localiser.dcm');
Concat_Localiser = [flip(Localiser(:,:,1,1),2);flip(Localiser(:,:,1,2),2);flip(Localiser(:,:,1,3),2)];
figure('color','w');
imagesc(Concat_Localiser,[0 0.6*max(Localiser,[],'all')])
colormap('gray'); axis image off

load('SatTFL_Reconstructed_B1Maps_Images.mat');
Concat_SatTFL = [squeeze(FAperV_Combined(:,:,16));squeeze(FAperV_Combined(:,16,:));flip(squeeze(FAperV_Combined(18,:,:)),1)];
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
Concat_SA2RAGE = [squeeze(FAperV_Combined(:,:,16));squeeze(FAperV_Combined(:,16,:));flip(squeeze(FAperV_Combined(18,:,:)),1)];
load('Sandwich_Reconstructed_B1Maps_Images.mat');
Concat_Sandwich = [squeeze(FAperV_Combined(:,:,16));squeeze(FAperV_Combined(:,16,:));flip(squeeze(FAperV_Combined(18,:,:)),1)];

Combined_Images = cat(2,Concat_SatTFL,Concat_SA2RAGE,Concat_Sandwich);

max_val = max([Concat_SatTFL,Concat_SA2RAGE,Concat_Sandwich],[],'all');

figure('color','w'); 
imagesc(Combined_Images,[0 max_val]); axis image off
colormap(turbo)
colorbar()


% Plot Lookup Table
figure('color','w','Name','Lookup Tables'); 
load('SA2RAGE/sa2rage_lookup_table.mat');
plot(x_query.*(180/pi),real(cos(x_query)),'k'); hold on
plot(x_query.*(180/pi),fx_interp,'r');
load('Sandwich/sandwich_lookup_table.mat');
plot(x_query.*(180/pi),fx_interp,'b'); hold on;
%load('SatTFL/sattfl_lookup_table.mat');
%plot(x_query,fx_interp,'k--'); hold on;

ylabel('Image Ratio, [a.u.]'); xlabel(['Flip Angle, [',char(176),']'])
ylim([-1 1]); xlim([0 180]);
legend('\it{Cosine}','SA2RAGE','Sandwich')

%% Plot Figure 2 Images

Image_n = 1:2;
Voltage_n = 4;
load('SatTFL_Reconstructed_B1Maps_Images.mat');
Concat_SatTFL = abs([squeeze(Mask(:,:,16).*Images(:,:,16,Image_n,Voltage_n));squeeze(Mask(:,16,:).*Images(:,16,:,Image_n,Voltage_n));flip(squeeze(Mask(18,:,:).*Images(18,:,:,Image_n,Voltage_n)),1)]);
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
Concat_SA2RAGE = abs([squeeze(Mask(:,:,16).*Images(:,:,16,Image_n,Voltage_n));squeeze(Mask(:,16,:).*Images(:,16,:,Image_n,Voltage_n));flip(squeeze(Mask(18,:,:).*Images(18,:,:,Image_n,Voltage_n)),1)]);
load('Sandwich_Reconstructed_B1Maps_Images.mat');
Concat_Sandwich = abs([squeeze(Mask(:,:,16).*Images(:,:,16,Image_n,Voltage_n));squeeze(Mask(:,16,:).*Images(:,16,:,Image_n,Voltage_n));flip(squeeze(Mask(18,:,:).*Images(18,:,:,Image_n,Voltage_n)),1)]);

Combined_Images = cat(2,reshape(Concat_SatTFL,[96 64]),reshape(Concat_SA2RAGE,[96 64]),reshape(Concat_Sandwich,[96 64]));

max_val = max(Combined_Images,[],'all');

figure('color','w'); 
imagesc(Combined_Images,[0 max_val]); axis image off
colorbar()
colormap(gray)
%% Plot Figure 3 Multiple Voltage Bland-Altman

cmap = turbo(9);
fig = figure('color','w','Units','Normalized','Position',[0,0,0.308333333333333,0.857407407407407]); tiledlayout(3,2,'padding','none','tilespacing','compact');
for Scheme_n = 1:3
if Scheme_n == 1
load('SatTFL_Reconstructed_B1Maps_Images.mat');
title_str1 = 'a) SatTFL'; title_str2 = 'b) SatTFL';
elseif Scheme_n == 2
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
title_str1 = 'c) SA2RAGE'; title_str2 = 'd) SA2RAGE';
elseif Scheme_n == 3
load('Sandwich_Reconstructed_B1Maps_Images.mat');
title_str1 = 'e) Sandwich'; title_str2 = 'f) Sandwich';
end
nexttile();    
for Voltage_n = size(Voltages,2):-1:1 % plot back to front, so we don't hide lower voltage measurements
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
if Scheme_n == 1
lgd = legend(cellstr(num2str(flip(Voltages,2)')),'Location','Northwest','NumColumns',3,'fontsize',6);
lgd.Title.String = 'Transmit Voltage [V]';
lgd.ItemTokenSize(1) = 3;
end
plot([0,180],[0,180],'k','HandleVisibility','off')
xlim([0 180]); ylim([0 180]);
xticks([0:20:180]); yticks([0:20:180]); 
axis square
xlabel(['Ground Truth Flip Angle, [',char(176),']'],'fontsize',8)
ylabel(['Measured Flip Angle, [',char(176),']'],'fontsize',8)
title(title_str1,'fontsize',10);
grid on
nexttile;
for Voltage_n = size(Voltages,2):-1:1
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[])-reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
plot([0,180],[0,0],'k','HandleVisibility','off')
xlim([0 180]); ylim([-20 20]);
xticks([0:20:180]); yticks([-20,0,20]);
axis square
xlabel(['Ground Truth Flip Angle, [',char(176),']'],'fontsize',8)
ylabel(['Measured - Ground Truth Flip Angle, [',char(176),']'],'fontsize',8)
title(title_str2,'fontsize',10);
grid on
end

%% Plot Figure 4 Single Voltage Bland-Altman

cmap = turbo(9);
Voltage_n = 4;
textx = 5; texty = 20;
PlotFontSize = 8;

fig = figure('color','w','Units','Normalized','Position',[0,0,0.308333333333333,0.857407407407407]); tiledlayout(3,2,'padding','none','tilespacing','compact');
for Scheme_n = 1:3
if Scheme_n == 1
load('SatTFL_Reconstructed_B1Maps_Images.mat');
title_str1 = 'a) SatTFL'; title_str2 = 'b) SatTFL';
elseif Scheme_n == 2
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
title_str1 = 'c) SA2RAGE'; title_str2 = 'd) SA2RAGE';
elseif Scheme_n == 3
load('Sandwich_Reconstructed_B1Maps_Images.mat');
title_str1 = 'e) Sandwich'; title_str2 = 'f) Sandwich';
end
nexttile();    
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
plot([0,180],[0,180],'k','HandleVisibility','off')
xlim([0 180]); ylim([0 180]);
xticks([0:20:180]); yticks([0:20:180]); 
axis square
xlabel(['Ground Truth Flip Angle, [',char(176),']'],'fontsize',PlotFontSize)
ylabel(['Measured Flip Angle, [',char(176),']'],'fontsize',PlotFontSize)
title(title_str1,'fontsize',PlotFontSize+2);
grid on
text(textx,180-9,['n = ',num2str(size(nonzeros(reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[])),1))],'Fontsize',PlotFontSize);


nexttile;
plot(reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[])-reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
plot([0,180],[0,0],'k','HandleVisibility','off')
xlim([0 180]); ylim([-20 20]);
xticks([0:20:180]); yticks([-20,0,20]);
axis square
xlabel(['Ground Truth Flip Angle, [',char(176),']'],'fontsize',PlotFontSize)
ylabel(['Measured - Ground Truth Flip Angle, [',char(176),']'],'fontsize',PlotFontSize)
title(title_str2,'fontsize',PlotFontSize+2);
grid on

Plot_mean = mean(reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]) - reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),'omitnan');
Plot_std = std(reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[]) - reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]),[],'all');
yline(Plot_mean,'--',num2str(Plot_mean,2),'LabelVerticalAlignment','middle','Fontsize',PlotFontSize); 
yline(Plot_mean + 1.96*Plot_std,':',num2str(Plot_mean + 1.96*Plot_std,2),'LabelVerticalAlignment','middle','Fontsize',PlotFontSize);
yline(Plot_mean - 1.96*Plot_std,':',num2str(Plot_mean - 1.96*Plot_std,2),'LabelVerticalAlignment','middle','Fontsize',PlotFontSize);
text(textx,texty-2,['RPC: ',num2str(1.96*Plot_std,2),char(176)],'Fontsize',PlotFontSize);
text(textx,texty-4,['CV: ',num2str(100*Plot_std/mean((reshape(Voltages(Voltage_n).*FAperV(:,:,:,Voltage_n),1,[])+reshape(Voltages(Voltage_n).*FAperV_Combined,1,[]))/2),2),'%'],'Fontsize',PlotFontSize);

end

%% Plot Figure 
cmap = turbo(9);
fig = figure('color','w','Units','Normalized','Position',[0.245833333333333,0.275,0.408854166666667,0.392592592592593]); tiledlayout(1,2,'padding','none','tilespacing','compact');

load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
SA2RAGE_FAperV = FAperV;
load('Sandwich_Reconstructed_B1Maps_Images.mat');
Sandwich_FAperV = FAperV;

nexttile();    
for Voltage_n = size(Voltages,2):-1:1 % plot back to front, so we don't hide lower voltage measurements
plot(reshape(Voltages(Voltage_n).*SA2RAGE_FAperV(:,:,:,Voltage_n),1,[]),reshape(Voltages(Voltage_n).*Sandwich_FAperV(:,:,:,Voltage_n),1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
lgd = legend(cellstr(num2str(flip(Voltages,2)')),'Location','Northwest','NumColumns',3,'fontsize',6);
lgd.Title.String = 'Transmit Voltage [V]';
lgd.ItemTokenSize(1) = 3;
plot([0,180],[0,180],'k','HandleVisibility','off')
xlim([0 180]); ylim([0 180]);
xticks([0:20:180]); yticks([0:20:180]); 
axis square
xlabel(['SA2RAGE Flip Angle, [',char(176),']'],'fontsize',8)
ylabel(['Sandwich Flip Angle, [',char(176),']'],'fontsize',8)
title('a)','fontsize',10);
grid on
nexttile;
for Voltage_n = size(Voltages,2):-1:1
plot(reshape(Voltages(Voltage_n).*0.5*(SA2RAGE_FAperV(:,:,:,Voltage_n)+Sandwich_FAperV(:,:,:,Voltage_n)),1,[]),reshape(Voltages(Voltage_n).*(SA2RAGE_FAperV(:,:,:,Voltage_n)-Sandwich_FAperV(:,:,:,Voltage_n)),1,[]),'.','color',cmap(Voltage_n,:),'Markersize',3,'HandleVisibility','off'); hold on
plot(-10,-10,'.','color',cmap(Voltage_n,:),'Markersize',20); hold on % legend tokens
end
plot([0,180],[0,0],'k','HandleVisibility','off')
xlim([0 180]); ylim([-20 20]);
xticks([0:20:180]); yticks([-20,0,20]);
axis square
xlabel(['Mean(SA2RAGE,Sandwich), [',char(176),']'],'fontsize',8)
ylabel(['SA2RAGE - Sandwich, [',char(176),']'],'fontsize',8)
title('b)','fontsize',10);
grid on

%% Mean FA per Volt across the entire head
load('SatTFL_Reconstructed_B1Maps_Images.mat');
disp(['SatTFL: Mean FA per Volt in whole head: ', num2str(mean(nonzeros(FAperV_Combined),'all')), ' ', char(177), ' ', num2str(mean(nonzeros(FAperV_Combined_SD),'all'))]);
disp(['SatTFL: Maximum FA per Volt in whole head: ', num2str(max(nonzeros(FAperV_Combined),[],'all')), '. Minimum ', num2str(min(nonzeros(FAperV_Combined),[],'all'))]);
disp(['SatTFL: Dynamic Range: ', num2str(max(nonzeros(FAperV_Combined),[],'all')./min(nonzeros(FAperV_Combined),[],'all'))]);

load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
disp(['SA2RAGE: Mean FA per Volt in whole head: ', num2str(mean(nonzeros(FAperV_Combined),'all')), ' ', char(177), ' ', num2str(mean(nonzeros(FAperV_Combined_SD),'all'))]);
disp(['SA2RAGE: Maximum FA per Volt in whole head: ', num2str(max(nonzeros(FAperV_Combined),[],'all')), '. Minimum ', num2str(min(nonzeros(FAperV_Combined),[],'all'))]);
disp(['SA2RAGE: Dynamic Range: ', num2str(max(nonzeros(FAperV_Combined),[],'all')./min(nonzeros(FAperV_Combined),[],'all'))]);

load('Sandwich_Reconstructed_B1Maps_Images.mat');
disp(['Sandwich: Mean FA per Volt in whole head: ', num2str(mean(nonzeros(FAperV_Combined),'all')), ' ', char(177), ' ', num2str(mean(nonzeros(FAperV_Combined_SD),'all'))])
disp(['Sandwich: Maximum FA per Volt in whole head: ', num2str(max(nonzeros(FAperV_Combined),[],'all')), '. Minimum ', num2str(min(nonzeros(FAperV_Combined),[],'all'))]);
disp(['Sandwich: Dynamic Range: ', num2str(max(nonzeros(FAperV_Combined),[],'all')./min(nonzeros(FAperV_Combined),[],'all'))]);




%% Relative Signal
Voltage_n = 4;
load('SA2RAGE_Reconstructed_B1Maps_Images.mat');
SA2RAGE_Images = abs(Mask.*Images(:,:,:,:,Voltage_n));
load('Sandwich_Reconstructed_B1Maps_Images.mat');
Sandwich_Images = abs(Mask.*Images(:,:,:,:,Voltage_n));

Ratio = SA2RAGE_Images./Sandwich_Images;

Ratio(isinf(Ratio)) = NaN;

mean(Ratio,1:3,'omitnan')

