%% Generate Masks
function Mask = Generate_Mask(Images)

% Normalise the reference images by maximum value at each transmit voltage
Ref_Images = Images(:,:,:,1,:)./max(abs(Images(:,:,:,1,:)),[],1:3);

% Binarise based on threshold value
threshold1 = 0.05;
Mask = squeeze(sum(abs(Ref_Images),5))./size(Ref_Images,5);
Mask(Mask < threshold1) = 0;
Mask(Mask > threshold1) = 1;

% Dilate and Erode mask
se = strel('square',3);
Mask = imdilate(Mask, se);
Mask = imerode(Mask, se);

Mask = logical(Mask); % Boolean

% Check Mask
%figure('color','w','Name','Images (Unmasked versus Masked)'); tiledlayout('flow');
%nexttile; imagesc(imtile(abs(squeeze(Images(:,:,:,1,1)))))
%nexttile; imagesc(imtile(abs(squeeze(permute(Mask,[1,2,3,5,4]).*Images(:,:,:,1,1)))))
end
