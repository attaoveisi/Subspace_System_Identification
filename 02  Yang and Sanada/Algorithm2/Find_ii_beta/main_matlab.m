clc
clear
close all
%%
for ii = 13:70
    for beta = 13:70
        
        [m,l,freq_vec,Gm,M,sysc,fmax] = Yang_Sanada_algorithm2(ii,beta);

        figure;
        k = 0;
        for i = 1:m
            for j = 1:l
                k = k+1;
                subplot(m,l,k)
                plot(freq_vec,20*log10(abs(squeeze(Gm(j,i,1:M)))));
                hold on
                [magG,phaseG,w] = bode(sysc(j,i),freq_vec*2*pi);
                for l = 1:size(w)
                    bb(l,1) = squeeze(magG(:,:,l));
                end
                plot(w/2/pi,20*log10(bb),'r')
                title(['ii = ',num2str(ii),' and beta = ',num2str(beta)])
                xlim([0 fmax])
            end
        end
    end
end