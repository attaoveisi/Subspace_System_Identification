clc
clear

m=1;                 % number of input
p=1;                 % number of output  

kk = 1;
for n = 16:16
    for q = 38:38
        [select_w,ff,Gjw,sysc] = q_n_finder(q,n,m,p);
        
        figure(kk);
        kk = kk+1;
        k = 0;
        for i = 1:m
            for j = 1:p
                k = k+1;
                subplot(m,p,k)
                plot(ff(select_w:end),20*log10(abs(squeeze(Gjw(j,i,select_w:end)))));
                title(['n is ',num2str(n),' and q is ',num2str(q)])
                hold on
                [magG,~]=freqresp(sysc(j,i),ff(select_w:end),'Hz');
                plot(ff(select_w:end),20*log10(abs(squeeze(magG))),'r')
                xlim([0 800])
                hold off
            end
        end
    end
end
