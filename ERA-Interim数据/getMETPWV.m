function [ pwv_met ] = getMETPWV( ZTD )
    [yp,yt]=readfile2;
    yt=yt+273.15;
    
    %TJCHվ��Ϣ
    B=31+17/60+5.84123/3600;
    L=121+29/60+51.84707/3600;
    %H=36.7305;
    H=3.54048420;

    %��Ϊ����
    B=B*pi/180;
    L=L*pi/180;
    
    ZHD=2.2779*yp/(1-0.00266*cos(2*B)-0.00028*H);
    Bigpi=10^5./(174073600./(0.72*(yt-273.15)+266.868)+7597.28);

    %��doy33 2�㵽doy45 14��
    ZHD=ZHD(3:304);
    Bigpi=Bigpi(3:304);
    pwv_met=(ZTD-ZHD').*Bigpi';

end

