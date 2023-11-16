function kdata = RORecenter(kdata,twix)

full_sz(1) = twix.image.centerLin(1)*2 -  2;
full_sz(2) = twix.image.centerPar(1)*2 -  2;
    
kc = twix.image.centerCol(1);
sz = size(kdata);
intended_ro = (sz(1) - kc/2)*2+1;
ro_add = intended_ro - sz(1);
kdata(intended_ro,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
kdata = circshift(kdata,[ro_add 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);

if(sz(3)<full_sz(1))
    kdata(:,:,full_sz(1),:,:,:,:,:,:,:,:,:,:,:,:,:) = 0;
end
if(sz(4)<full_sz(2))
    kdata(:,:,:,full_sz(2),:,:,:,:,:,:,:,:,:,:,:,:) = 0;
end