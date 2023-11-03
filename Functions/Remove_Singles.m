%% Functions
function [Exclusion_Mask] = Remove_Singles(Exclusion_Mask)
% Remove values which are:
% (a) part of a very small island, and 
% (b) duplicate (i.e. are included from >1 dataset). 
max_sz = 10;
BW = false(size(Exclusion_Mask));
for Voltage_n = 1:size(Exclusion_Mask,4)
duplicate_locs = sum(Exclusion_Mask,4) < 2;
BW(:,:,:,Voltage_n) = ~bwareaopen(~Exclusion_Mask(:,:,:,Voltage_n),max_sz,6); 
Temp = Exclusion_Mask(:,:,:,Voltage_n);
Temp(logical((BW(:,:,:,Voltage_n)-Exclusion_Mask(:,:,:,Voltage_n)).*duplicate_locs)) = 1;
Exclusion_Mask(:,:,:,Voltage_n) = Temp; clear Temp
end
end