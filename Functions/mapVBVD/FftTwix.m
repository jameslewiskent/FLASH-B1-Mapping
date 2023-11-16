function [im] = FftTwix(kdata)
    im = fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(kdata,1),3),4),[],1),[],3),[],4),1),3),4);
end

