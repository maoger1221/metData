function [press_m,temp_m]=readfile2
    data=load('met_all.txt');
    k=0;sump=0;n=1;sumt=0;
    for i=2:length(data)
        if data(i,4)==data(i-1,4)
            sump=sump+data(i,7);
            sumt=sumt+data(i,8);
            k=k+1;
        else
            press_m(n)=sump/k;
            temp_m(n)=sumt/k;
            n=n+1;
            k=1;
            sump=data(i,7);
            sumt=data(i,8);
        end
    end
end