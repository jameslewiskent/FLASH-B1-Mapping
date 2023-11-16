%% Reconstruct Data from TWIX
% Please contact james.kent@ndcn.ox.ac.uk for access to raw Twix data
clearvars
close all
main_directory = fileparts(matlab.desktop.editor.getActiveFilename);
cd(main_directory); addpath(genpath(main_directory));
n = 2; kernel = [5,5]; eig_thresh = 0.02;
Desired_ro = 32; Desired_pe = 32;
redo_image_recon = 1; % Re-run the image reconstruction
%% 2D MS SatTFL Recon
clearvars -except n main_directory kernel eig_thresh Desired_ro Desired_pe Voltages redo_image_recon
cd([main_directory,filesep,'SatTFL',filesep,'20231101b'])
MIDs = [84,85,86,87,88,89,90,91,92];
if ~isfile('SatTFL_Reconstructed_B1Maps_Images.mat') || redo_image_recon == 1
    Images = zeros(Desired_ro,Desired_pe,32,2,length(MIDs));
    parfor MID_n = 1:length(MIDs)
        MID = MIDs(MID_n);
        twix = mapVBVD(MID);
        twix{1,n}.image.flagRemoveOS = true;
        Voltages(MID_n) = twix{1,n}.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;

        kdata = squeeze(twix{1,n}.image(''));
        dims = twix{1,n}.image.sqzDims;
        
        % Sort slices
        [~, C] = sort(twix{1,n}.hdr.Config.chronSliceIndices(1:twix{1,n}.image.NSli));
        kdata = kdata(:,:,:,C,:);
        for Slice_n = 1:32
            kdata_slice = squeeze(kdata(:,:,:,Slice_n,:));
            rx_sens = tx_espirit(permute(kdata_slice,[1,3,4,2]),[Desired_ro,Desired_pe], kernel, eig_thresh);
            
            if size(kdata_slice,1) < Desired_ro || size(kdata_slice,2) < Desired_pe
                % Zero-fill and Hanning filter
                kdata_zp = squeeze(Zero_Fill(twix{1,n},kdata_slice,Desired_ro,Desired_pe));
            else
                kdata_zp = kdata_slice;
            end
            
            % fft in col and lin dimension:
            fft_dims = [1 3];
            for f = fft_dims
                kdata_zp = ifftshift(ifft(fftshift(kdata_zp,f),[],f),f);
            end
            
            % Combine sensitivity maps with images, sum across receive channels
            Images(:,:,Slice_n,:,MID_n) = squeeze( sum( permute(kdata_zp,[1 3 2 4]).*conj(rx_sens) ,3) );
        end
    end
    % Check Images
    figure('color','w','Name','SatTFL Images'); tiledlayout('flow');
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,1,:))))); axis image off
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,2,:))))); axis image off
    
    % Flip SatTFL images to be in same orientation as other datasets
    Images = flip(Images,3);    
    
    save('SatTFL_Reconstructed_B1Maps_Images.mat','Images','Voltages');
else
    load('SatTFL_Reconstructed_B1Maps_Images.mat','Images','Voltages');
end

Mask = Generate_Mask(Images);
save('SatTFL_Reconstructed_B1Maps_Images.mat','Mask','-append');

Maps = zeros(Desired_ro,Desired_pe,32,length(MIDs));
for MID_n = 1:length(MIDs)
    Maps(:,:,:,MID_n) = 180/pi.*acos(real(Images(:,:,:,2,MID_n)./Images(:,:,:,1,MID_n)));
end
save('SatTFL_Reconstructed_B1Maps_Images.mat','Maps','-append');

% Check B1 Maps
figure('color','w','Name','SatTFL B1 Maps'); imagesc(imtile(abs(squeeze(Maps(:,:,16,:)))),[0 150]); axis image off

% Combine multi-voltage B1 maps
[FAperV,FAperV_Combined,FAperV_Combined_SD] = Combine_Multivoltage_Datasets(Images,Maps,Mask,Voltages,1:5); % Do not include clipped saturation pulses (> 275 V)
save('SatTFL_Reconstructed_B1Maps_Images.mat','FAperV','FAperV_Combined','FAperV_Combined_SD','-append');


%% 3D SA2RAGE Recon
clearvars -except n main_directory kernel eig_thresh Desired_ro Desired_pe Voltages redo_image_recon
cd([main_directory,filesep,'SA2RAGE',filesep,'20231101b']);
MIDs = [93,94,95,96,97,98,99,100,101];

if ~isfile('SA2RAGE_Reconstructed_B1Maps_Images.mat') || redo_image_recon == 1
    Images = zeros(Desired_ro,Desired_pe,32,2,length(MIDs));
    parfor MID_n = 1:length(MIDs)
        MID = MIDs(MID_n);
        twix = mapVBVD(MID);
        twix{1,n}.image.flagRemoveOS = true;
        Voltages(MID_n) = twix{1,n}.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;
        
        kdata = squeeze(twix{1,n}.image(''));
        dims = twix{1,n}.image.sqzDims;
        
        kdata = flip(kdata,5); % Flip SA2RAGE to put 'reference' image first for convenience
        
        % FFT in par dimension
        for f = 4
            kdata = ifftshift(ifft(fftshift(kdata,f),[],f),f);
        end
        
        for Part_n = 1:32
            kdata_slice = squeeze(kdata(:,:,:,Part_n,:));
            rx_sens = tx_espirit(permute(kdata_slice,[1,3,4,2]),[Desired_ro,Desired_pe], kernel, eig_thresh);
            
            if size(kdata_slice,1) < Desired_ro || size(kdata_slice,2) < Desired_pe
                kdata_zp = squeeze(Zero_Fill(twix{1,n},kdata_slice,Desired_ro,Desired_pe)); % Zero-fill and Hanning filter
            else
                kdata_zp = kdata_slice;
            end
            
            % fft in col and lin dimension:
            fft_dims = [1 3];
            for f = fft_dims
                kdata_zp = ifftshift(ifft(fftshift(kdata_zp,f),[],f),f);
            end
            
            % Combine sensitivity maps with images, sum across receive channels
            Images(:,:,Part_n,:,MID_n) = squeeze( sum( permute(kdata_zp,[1 3 2 4]).*conj(rx_sens) ,3) );
        end
    end
    % Check Images
    figure('color','w','Name','SA2RAGE Images'); tiledlayout('flow');
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,1,:))))); axis image off
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,2,:))))); axis image off
    
    save('SA2RAGE_Reconstructed_B1Maps_Images.mat','Images','Voltages');
else
    load('SA2RAGE_Reconstructed_B1Maps_Images.mat','Images','Voltages');
end

Mask = Generate_Mask(Images);
save('SA2RAGE_Reconstructed_B1Maps_Images.mat','Mask','-append');

Maps = zeros(Desired_ro,Desired_pe,32,length(MIDs));
load('../sa2rage_lookup_table.mat'); % Load lookup table
for MID_n = 1:length(MIDs)
    Image_Ratio = real(Images(:,:,:,2,MID_n)./Images(:,:,:,1,MID_n));
    Image_Ratio(Image_Ratio > max(fx_interp)) = max(fx_interp);
    Image_Ratio(Image_Ratio < min(fx_interp)) = min(fx_interp);
    for x = 1:size(Image_Ratio,1)
        for y = 1:size(Image_Ratio,2)
            for z = 1:size(Image_Ratio,3)
                [~,min_ind] = min(abs(Image_Ratio(x,y,z) - fx_interp));
                Maps(x,y,z,MID_n) = (180/pi)*x_query(min_ind); % Measured FA in degrees
            end
        end
    end
end
save('SA2RAGE_Reconstructed_B1Maps_Images.mat','Maps','-append');

% Check B1 Maps
figure('color','w','Name','SA2RAGEs B1 Maps'); imagesc(imtile(abs(squeeze(Maps(:,:,16,:)))),[0 150]); axis image off

% Combine multi-voltage B1 maps
[FAperV,FAperV_Combined,FAperV_Combined_SD] = Combine_Multivoltage_Datasets(Images,Maps,Mask,Voltages,1:9);
save('SA2RAGE_Reconstructed_B1Maps_Images.mat','FAperV','FAperV_Combined','FAperV_Combined_SD','-append');

%% 3D Sandwich Recon
clearvars -except n main_directory kernel eig_thresh Desired_ro Desired_pe Voltages redo_image_recon
cd([main_directory,filesep,'Sandwich',filesep,'20231101b']);
MIDs = [102,103,104,105,106,107,108,109,110];
if ~isfile('Sandwich_Reconstructed_B1Maps_Images.mat') || redo_image_recon == 1
    Images = zeros(Desired_ro,Desired_pe,32,2,length(MIDs));
    parfor MID_n = 1:length(MIDs)
        MID = MIDs(MID_n);
        twix = mapVBVD(MID);
        twix{1,n}.image.flagRemoveOS = true;
        Voltages(MID_n) = twix{1,n}.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;
        
        kdata = squeeze(twix{1,n}.image(''));
        dims = twix{1,n}.image.sqzDims;
        
        % FFT in par dimension
        for f = 4
            kdata = ifftshift(ifft(fftshift(kdata,f),[],f),f);
        end
        
        for Part_n = 1:32
            kdata_slice = squeeze(kdata(:,:,:,Part_n,:));
            rx_sens = tx_espirit(permute(kdata_slice,[1,3,4,2]),[Desired_ro,Desired_pe], kernel, eig_thresh);
            
            if size(kdata_slice,1) < Desired_ro || size(kdata_slice,2) < Desired_pe
                kdata_zp = squeeze(Zero_Fill(twix{1,n},kdata_slice,Desired_ro,Desired_pe)); % Zero-fill and Hanning filter
            else
                kdata_zp = kdata_slice;
            end
            
            % fft in col and lin dimension:
            fft_dims = [1 3];
            for f = fft_dims
                kdata_zp = ifftshift(ifft(fftshift(kdata_zp,f),[],f),f);
            end
            
            % Combine sensitivity maps with images, sum across receive channels
            Images(:,:,Part_n,:,MID_n) = squeeze( sum( permute(kdata_zp,[1 3 2 4]).*conj(rx_sens) ,3) );
        end
    end
    % Check Images
    figure('color','w','Name','Sandwich Images'); tiledlayout('flow');
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,1,:))))); axis image off
    nexttile; imagesc(imtile(abs(squeeze(Images(:,:,16,2,:))))); axis image off
    
    save('Sandwich_Reconstructed_B1Maps_Images.mat','Images','Voltages');
else
    load('Sandwich_Reconstructed_B1Maps_Images.mat','Images','Voltages');
end

Mask = Generate_Mask(Images);
save('Sandwich_Reconstructed_B1Maps_Images.mat','Mask','-append');

Maps = zeros(Desired_ro,Desired_pe,32,length(MIDs));
load('Sandwich/sandwich_lookup_table.mat'); % Load lookup table
for MID_n = 1:length(MIDs)
    Image_Ratio = real(Images(:,:,:,2,MID_n)./Images(:,:,:,1,MID_n));
    Image_Ratio(Image_Ratio > max(fx_interp)) = max(fx_interp);
    Image_Ratio(Image_Ratio < min(fx_interp)) = min(fx_interp);
    for x = 1:size(Image_Ratio,1)
        for y = 1:size(Image_Ratio,2)
            for z = 1:size(Image_Ratio,3)
                [~,min_ind] = min(abs(Image_Ratio(x,y,z) - fx_interp));
                Maps(x,y,z,MID_n) = (180/pi)*x_query(min_ind); % Measured FA in degrees
            end
        end
    end
end
save('Sandwich_Reconstructed_B1Maps_Images.mat','Maps','-append');

% Check B1 Maps
figure('color','w','Name','Sandwich B1 Maps'); imagesc(imtile(abs(squeeze(Maps(:,:,16,:)))),[0 150]); axis image off

% Combine multi-voltage B1 maps
[FAperV,FAperV_Combined,FAperV_Combined_SD] = Combine_Multivoltage_Datasets(Images,Maps,Mask,Voltages,1:9);
save('Sandwich_Reconstructed_B1Maps_Images.mat','FAperV','FAperV_Combined','FAperV_Combined_SD','-append');
