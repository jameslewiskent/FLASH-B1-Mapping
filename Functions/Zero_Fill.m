function [kdata] = Zero_Fill(twix,kdata,Desired_ro,Desired_pe,Desired_pa)
% Zero-fill k-space data in read out and phase encode
sz = size(twix.image('')); % pre-filled size
if twix.image.flagRemoveOS == true
    os = 1;
else
    os  = 2;
end

intended_ro =  ((sz(1)/os - twix.image.centerCol(1)/2)*2+ 1);
if nargin <3  || Desired_ro < intended_ro
Desired_ro = intended_ro;
end 
if Desired_ro ~= sz(1)
ro_add = intended_ro - sz(1) + os.*(Desired_ro-intended_ro)/2;
kdata = kdata.*hann(size(kdata,1)); % Hanning filter
kdata(Desired_ro,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[ro_add 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
end


intended_pe = 2*(sz(3) - twix.image.centerLin(1) + 1);
if nargin <4 || Desired_pe < intended_pe
    Desired_pe = intended_pe;
end
if Desired_pe ~= sz(3)
pe_add = intended_pe - sz(3) + (Desired_pe-intended_pe)/2;
kdata = kdata.*permute(hann(size(kdata,3)),[2 3 1]); % Hanning filter
kdata(:,:,Desired_pe,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[0 0 pe_add 0 0 0 0 0 0 0 0 0 0 0 0 0]);
end


intended_pa = 2*(sz(4) - twix.image.centerPar(1) + 1);
if nargin <5 || Desired_pa < intended_pa
    Desired_pa = intended_pa;
end
if Desired_pa ~= sz(4) && sz(4) ~= 1 % if not the desired partition size and the size of the partition isn't currently 1 (otherwise doubles to 2)
pa_add = intended_pa - sz(4) + (Desired_pa-intended_pa)/2;
kdata = kdata.*permute(hann(size(kdata,4)),[2 3 4 1]); % Hanning filter
kdata(:,:,:,Desired_pa,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[0 0 0 pa_add 0 0 0 0 0 0 0 0 0 0 0 0]);
end

end