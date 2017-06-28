function [ pwv_ncep,yp,yt ] = getNCEP(ZTD)%����ֵʱ�䷶Χdoy33 0�㵽doy45 24�㣬��53��ʱ�䣨13��*һ��4��+1��doy46 0�㣩��
%(30,120) 1 (30,122.5) 2
%(32.5,120) 3 (32.5,122.5) 4
%120E122.5E����Ӧ���к�49��50
%30N32.5N����Ӧ���к�25��24
tid=netcdf.open('air.sig995.2015.nc','NOWRITE');
pid=netcdf.open('pres.sfc.2015.nc','NOWRITE');
hid=netcdf.open('hgt.sfc.nc','NOWRITE');

%datestr(datenum(1800,1,1)+����)
%ncep��ʱ���1800.1.1��ʼ(ʱ���ΪUTC)
time = netcdf.getVar(tid,2,'double');
%2015.2.2��doy33����1800.1.1��Сʱ��
hours1=(datenum(2015,2,2)-datenum(1800,1,1))*24;
%2015.2.15��doy46����1800.1.1��Сʱ��
hours2=(datenum(2015,2,15)-datenum(1800,1,1))*24;
%����2015.2.2 0���Ӧ���к�
rows1=find(time==hours1);
%����2015.2.15 0���Ӧ���к�
rows2=find(time==hours2);
%��r����doy33 0�㵽doy45 24��
r=rows1:rows2;
%r������
size2=size(r,2);

air = netcdf.getVar(tid,3,'double');
t1=air(49,25,r);
t2=air(50,25,r);
t3=air(49,24,r);
t4=air(50,24,r);
%����ά����ά
temp=[reshape(t1,size2,1) reshape(t2,size2,1) reshape(t3,size2,1) reshape(t4,size2,1)];


pres = netcdf.getVar(pid,3,'double');
p1=pres(49,25,r);
p2=pres(50,25,r);
p3=pres(49,24,r);
p4=pres(50,24,r);
%����ά����ά
press=[reshape(p1,size2,1) reshape(p2,size2,1) reshape(p3,size2,1) reshape(p4,size2,1)];


hgt = netcdf.getVar(hid,3,'double');
h1=hgt(49,25);
h2=hgt(50,25);
h3=hgt(49,24);
h4=hgt(50,24);
%����ά����ά
height=[h1 h2 h3 h4];

%TJCHվ��Ϣ
B=31+17/60+5.84123/3600;
L=121+29/60+51.84707/3600;
%H=36.7305;
H=3.54048420;

%��ѹ��ֱ��ֵ
p1v=press(:,1).*exp(-(H-height(1))*9.8*0.0289/8.31./temp(:,1));
p2v=press(:,2).*exp(-(H-height(2))*9.8*0.0289/8.31./temp(:,2));
p3v=press(:,3).*exp(-(H-height(3))*9.8*0.0289/8.31./temp(:,3));
p4v=press(:,4).*exp(-(H-height(4))*9.8*0.0289/8.31./temp(:,4));
%�¶ȴ�ֱ��ֵ
t1v=temp(:,1)+0.98*(H-height(1))/100;
t2v=temp(:,2)+0.98*(H-height(2))/100;
t3v=temp(:,3)+0.98*(H-height(3))/100;
t4v=temp(:,4)+0.98*(H-height(4))/100;

%��Ϊ����
B=B*pi/180;
L=L*pi/180;
Bj=[30*pi/180 30*pi/180 32.5*pi/180 32.5*pi/180];
Lj=[120*pi/180 122.5*pi/180 120*pi/180 122.5*pi/180];

wj=(acos(sin(Bj)*sin(B)+cos(Bj)*cos(B).*cos(Lj-L))*6371004).^(-2);
w1=wj(1)/sum(wj);
w2=wj(2)/sum(wj);
w3=wj(3)/sum(wj);
w4=wj(4)/sum(wj);

p=p1v*w1+p2v*w2+p3v*w3+p4v*w4;
t=t1v*w1+t2v*w2+t3v*w3+t4v*w4;
p=p/100;%��λ����hPa
x1=1:size(p,1);

xx=1:1/6:size(p,1);%��ֵ�����Ϊ1Сʱ��
y1=p;
y2=t;
yp=spline(x1,y1,xx);
yt=spline(x1,y2,xx);

ZHD=2.2779*yp/(1-0.00266*cos(2*B)-0.00028*H);
Bigpi=10^5./(174073600./(0.72*(yt-273.15)+266.868)+7597.28);

%��doy33 2�㵽doy45 14��
ZHD=ZHD(3:304);
Bigpi=Bigpi(3:304);
pwv_ncep=(ZTD-ZHD').*Bigpi';
end

