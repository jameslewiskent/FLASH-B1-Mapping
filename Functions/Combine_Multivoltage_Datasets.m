function [FAperV,FAperV_Combined,FAperV_Combined_SD] = Combine_Multivoltage_Datasets(Images,Maps,Mask,Voltages,Voltage_Indices)
% Code adapted from Matthijs de Buck. to handle >3 voltage datasets
% Voltage_Indices dictates which voltages are used in the final combined FAperV

% Check voltage values are increasing
if any(diff(Voltages) <= 0)
    error('Maps are not ordered in increasing voltage!')
end

Reference_Images =  squeeze(abs(Images(:,:,:,1,:)));
FAperV =  abs(Mask.*Maps(:,:,:,:)./permute(Voltages,[1 3 4 2]));

figure('color','w','Name','Flip Angle per Volt'); imagesc(imtile(squeeze(FAperV(:,:,16,:))),[0 max(FAperV,[],'all')]); axis image off

Low_Sig_Mask = false([size(Reference_Images,1:3),(size(Reference_Images,4)-1)]); % Pre-allocate low signal mask
for Voltage_n = 1:(size(Reference_Images,4)-1)
    Low_Sig_Mask(:,:,:,Voltage_n) = (Reference_Images(:,:,:,Voltage_n)<0.15*max(Reference_Images(:,:,:,end))) .* (Reference_Images(:,:,:,Voltage_n)*1.5 < Reference_Images(:,:,:,Voltage_n + 1)) .* (FAperV(:,:,:,Voltage_n) > FAperV(:,:,:,Voltage_n + 1));
end
Reference_Images = Reference_Images./max(Reference_Images,[],1:3); % Normalise
SD = std(Reference_Images(~Mask));

% Criterion 1: There has to be a minimum amount of signal in the reference images
Exclusion_Mask = Reference_Images.*Mask < SD;
Exclusion_Mask(cat(4,Low_Sig_Mask,false(size(Low_Sig_Mask,1:3)))) = true; % Partial fix for low-intensity regions
Exclusion_Mask = Remove_Singles(Exclusion_Mask);

figure('Name','Criterion 1'); imagesc(imtile(squeeze(Exclusion_Mask(:,:,16,:)))); axis image off

% Criterion 2: Higher-Voltage FA map should have a flip angle greater than any of the lower-voltage maps
for Voltage_n = 2:size(Reference_Images,4)
    Temp = Exclusion_Mask(:,:,:,Voltage_n);
    for Lower_Voltage_n = 1:(Voltage_n-1)
    Temp(Voltages(Voltage_n).*(FAperV(:,:,:,Voltage_n)) <= Voltages(Lower_Voltage_n).*(FAperV(:,:,:,Lower_Voltage_n) .* ~Exclusion_Mask(:,:,:,Lower_Voltage_n))) = 1;
    end
    Exclusion_Mask(:,:,:,Voltage_n) = Temp; clear Temp
end
Exclusion_Mask = Remove_Singles(Exclusion_Mask);

figure('Name','Criterion 2'); imagesc(imtile(squeeze(Exclusion_Mask(:,:,16,:)))); axis image off

% Combine maps
FAperV_Masked = FAperV;
FAperV_Masked(Exclusion_Mask) = NaN;
FAperV_Combined = mean(FAperV_Masked(:,:,:,Voltage_Indices),4,'omitnan');
FAperV_Combined(isnan(FAperV_Combined)) = 0;
FAperV_Combined_SD = std(FAperV_Masked,[],4,'omitnan');
FAperV_Combined_SD(isnan(FAperV_Combined_SD)) = 0;

figure('color','w','Name','Combined Flip Angle per Volt'); tiledlayout('flow');
nexttile; imagesc(imtile(FAperV_Combined),[0 max(FAperV_Combined,[],'all')]); title('Mean'); axis image off;
nexttile; imagesc(imtile(FAperV_Combined_SD),[0 max(FAperV_Combined,[],'all')]); title('SD'); axis image off

end
