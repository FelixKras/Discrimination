ind=[0:2^16-1];
Red=[zeros(1,1000) linspace(0,255,25e3-1000) 255*ones(1,2^16-25e3)];
Green=[linspace(0,64,2^15) 64*ones(1,2^16-2^15)];
Blue=[linspace(0,128,700) zeros(1,2^15-700) 128*ones(1,2^16-2^15)];
dlmwrite('c:\paletteforGalit.txt',[uint16(ind); uint16(Red); uint16(Green); uint16(Blue)]','delimiter',' ');