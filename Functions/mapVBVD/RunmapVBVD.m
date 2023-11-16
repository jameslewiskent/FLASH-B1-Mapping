%addpath('C:\Users\jkent\OneDrive - Nexus365\Scanner Data\Raw Data\23032021')
% clc
% clear all

twix = mapVBVD(52);
disp(['TR = ', num2str( twix.hdr.Config.TR /1e6),' s'])
disp(['FlipAngle = ', num2str(twix.hdr.Meas.adFlipAngleDegree)])
%disp(['PrepFlipAngle = ', num2str(twix.hdr.Config.PrepFlipAngle)])
disp(['Tx Ref. Amp = ', num2str( twix.hdr.Meas.TransmitterReferenceAmplitude),' V'])
disp(['BW = ', num2str( 0),' Hz/pixel'])
disp(['NSeg = ', num2str(twix.image.NSeg)])
disp(['EPIFactor = ', num2str(twix.hdr.Meas.lEPIFactor)])
disp(['EnableRFSpoiling (1 is on) = ', num2str(twix.hdr.Meas.ulEnableRFSpoiling)])
disp(['PhaseResolution = ', num2str(twix.hdr.Meas.dPhaseResolution)])
disp(['PhasePartialFourierFactor = ', num2str(twix.hdr.Meas.PhasePartialFourierFactor)])

Prep_FA = twix.hdr.Meas.adFlipAngleDegree *10;
kdata = twix.image('');
size(kdata)

% [rf,in] = twix.vop.readVOP();
% data = rf{1}(:,:,:).*conj(rf{1}(:,1,:));
% figure()
% for n = 1:2:16
% plot(squeeze(abs(data(250,n,:))),'x-')
% %xlim([0 3])
% xlabel('Prep-Pulse #')
% ylabel('Magnitude')
% hold on
% end
% legend('Tx0','Tx1','Tx2','Tx3','Tx4','Tx5','Tx6','Tx7')
% 
% figure()
% for n = 1:2:16
% hold on
% plot(squeeze(angle(data(250,n,:)).*180/pi),'x-')
% %xlim([0 3])
% xlabel('RF-Pulse #')
% ylabel('Phase')
% end
%% View segments
twix.image.flagIgnoreSeg = false;
twix.image.flagRemoveOS = true;
kdata = twix.image('');

[kdata] = zero_filling(twix,kdata);

part = 1;
Ave = 1;
Phs = 1;
Eco  =1;
Rep = 1;
for freeloop = 1:twix.image.NSet
    figure('Name',['K-space Data',num2str(freeloop)],'Position',[556.2 250.6  552.8  528])
    for slice = 1
        for RxChannel = 1
            for Seg_n = 1:twix.image.NSeg
                subplot(1,twix.image.NSeg,Seg_n)
                %                                          Col RxCha Lin Par Sli Ave Phs Eco Rep Set Seg
                imagesc((abs(squeeze(kdata(:,RxChannel,:,part,slice,Ave,Phs,Eco,Rep,freeloop,Seg_n))).^0.2),[0 max(abs(kdata),[],'all')^0.2]);
                axis normal off
            end
        end
    end
end



%% View k-space of receive channels with segmentation ignored
twix.image.flagIgnoreSeg = true;
twix.image.flagRemoveOS =true;
kdata = twix.image('');

[kdata] = zero_filling(twix,kdata);

figure('Name','','Position',[556.2 250.6  552.8  528])
imagesc(imtile((abs(squeeze(permute(kdata(:,:,:,1,1,1,1,1,1,1,1),[1 3 2]))).^0.2),'Gridsize',[3,10]), [0 max(abs(kdata),[],'all')^0.2]);
axis normal off

figure('Name','','Position',[556.2 250.6  552.8  528])
imagesc(imtile((abs(squeeze(permute(kdata(:,:,:,1,1,1,1,1,1,2,1),[1 3 2]))).^0.2),'Gridsize',[3,10]), [0 0.4130]);
axis normal off

figure('Name','','Position',[556.2 250.6  552.8  528])
imagesc(imtile((angle(squeeze(permute(kdata(:,:,:,1,1,1,1,1,1,1,1),[1 3 2]))))));
axis normal off
figure('Name','','Position',[556.2 250.6  552.8  528])

imagesc(imtile((angle(squeeze(permute(kdata(:,:,:,1,1,1,1,1,1,2,1),[1 3 2]))))));
axis normal off

%% Construct Image
twix.image.flagIgnoreSeg = true;
twix.image.flagRemoveOS =true;

os = 2; % oversampling factor
nMID = 1;

% read in the data
kdata = squeeze(twix.image(:,:,:,1,1,1,1,1,1,:));
kdata = zero_filling(twix,kdata);

% Calculate sensitivity profiles from ESPIRIT
sens = tx_espirit(permute(kdata,[1 3 4 2]),[64 64]);
%
% fft in col and lin dimension:
fft_dims = [1 3];
for f = fft_dims
    kdata = ifftshift(ifft(fftshift(kdata,f),[],f),f);
end
sep_coil_imgs = kdata.*repmat(conj(permute(sens,[1 3 2])),[1 1 1 2]);

img= squeeze(sum(sep_coil_imgs,2));

B1Map = real(acos(((img(:,:,2)./img(:,:,1))))).*180/pi;

if ~exist('Mask','var')
    Mask = mask_data(img);
end
        %B1Map = Mask.*B1Map;
       %img =  Mask.*img;

figure()
clims = [min(min(abs(img(:,:,1)),[],'all'),min(abs(img(:,:,2)),[],'all')), max(max(abs(img(:,:,1)),[],'all'),max(abs(img(:,:,2)),[],'all'))];
%clims = [[0,0.00019954948]];
subplot_tight(2,2,1)
imagesc((abs(img(:,:,1))),clims)
title('Reference')
axis equal off
subplot_tight(2,2,2)
imagesc((abs(img(:,:,2))),clims)
title('Prep')
%colormap gray
axis equal off
subplot_tight(2,2,3.5)
imagesc(B1Map,[0 120]) %
title('B1 Map')
h= colorbar;
ylabel(h, 'Measured \alpha, degrees')
%colormap gray
axis equal off

%% View Individual Channel Images (+Sensitivities)

figure('Name','Rx Coil Data IT1')
imagesc(imtile((abs(squeeze(permute(sep_coil_imgs(:,:,:,1),[1 3 2 4]))))));
axis square off

figure('Name','Rx Coil Data IT2')
imagesc(imtile((abs(squeeze(permute(sep_coil_imgs(:,:,:,2),[1 3 2 4]))))),[10e-8 10e-6]);
axis square off

imagesc(abs(squeeze(sep_coil_imgs(:,19,:,2))),[5e-10 5e-6])


%% Read in VOP

[rf,in] = twix.vop.readVOP();

plot(squeeze(abs(rf{2}(250,1,:))))

%% Function Definitions
function [mask] = mask_data(img)

Reference_img = abs(img(:,:,1));
mask = zeros(size(Reference_img));

mask(Reference_img >2e-5) = 1; %prctile(Reference_img,64.45,'all')
%histogram(Reference_img,100)
end

function [kdata] = zero_filling(twix,kdata)
% zero-filling
%twix.image.flagIgnoreSeg = true;
% twix.image.flagRemoveOS = true;
%kdata = twix.image('');
sz = size(twix.image('')); % pre-filled size

kc = twix.image.centerCol(1);

%full_sz(1) = 50;%twix.image.centerLin(1)+1;
full_sz(2) = twix.image.centerPar(1)*2 - 2;


intended_ro = (sz(1) - kc/2)*2+1;
ro_add = intended_ro - sz(1);
kdata(intended_ro,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[ro_add 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);


intended_pe = 2*(sz(3) - twix.image.centerLin(1) + 1);
Full_pe = 64;
pe_add = intended_pe - sz(3) + (Full_pe-intended_pe)/2;
kdata(:,:,Full_pe,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[0 0 pe_add 0 0 0 0 0 0 0 0 0 0 0 0 0]);



%         if(sz(3)<full_sz(1))
%             kdata(:,:,full_sz(1),:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
%         end
%         if(sz(4)<full_sz(2))
%             kdata(:,:,:,full_sz(2),:,:,:,:,:,:,:,:,:,:,:,:) = 0;
%         end
end
