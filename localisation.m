function sortie_loca=localisation(gauss,signal,wav_name,cas_de_figure)

fs=44100;
sortie_space=spatialisation(gauss,signal,wav_name);
duree_son = length(signal) / fs;
outputSignal=sortie_space.outputSignal;
N=sortie_space.N;
D=sortie_space.D;
Lframe=sortie_space.Lframe;
Taille_de_1_position=sortie_space.Taille_de_1_position;
freqIndexes=sortie_space.freqIndexes;
%Taille_groupe=sortie_space.Taille_groupe;
B=sortie_space.B;
Ntheta=sortie_space.Ntheta;
P=sortie_space.P;
sigma=sortie_space.sigma;
thetaArg=sortie_space.thetaArg;
azimuth=sortie_space.azimuth;


x1=outputSignal(:,1);
x2=outputSignal(:,2);
% decoupage en N intervalles de taille Lframe
w=hanning(Lframe);
%w=1;
coef=sqrt(N/duree_son);
%%

for numero_exp=1:length(D)
    deb=Taille_de_1_position*(numero_exp-1)+round(Taille_de_1_position/2)-round(Lframe*2.5/2); %debut a peu pres au milieu dun sig localise
    % segments extremites
    %fprintf('deb: %d\n',deb)
    %fprintf('%d:%d\n',deb+1,deb+Lframe);
    x1t{1}(:,1)=x1(deb+1:deb+Lframe).*w;
    x2t{1}(:,1)=x2(deb+1:deb+Lframe).*w;
%     X1{1}=coef*fft(x1t{1});
%     X2{1}=coef*fft(x2t{1});
    X1{1}=coef*fft(x1t{1});
    X2{1}=coef*fft(x2t{1});
    
    
    for i=2:4 % segments internes
        %fprintf('%d:%d\n',deb+Lframe*(i-1)*.75+1,deb+Lframe*(i-1)*.75+Lframe);
        x1t{i}(:,1)=x1(deb+floor(Lframe/2)*(i-1)+1:deb+floor(Lframe/2)*(i-1)+Lframe).*w;
        x2t{i}(:,1)=x2(deb+floor(Lframe/2)*(i-1)+1:deb+floor(Lframe/2)*(i-1)+Lframe).*w;
        %         % fft
        X1{i}=coef*fft(x1t{i});
        X2{i}=coef*fft(x2t{i});
    end
    for k=1:B
        for i=1:4
            Zint{k}(2*i-1:2*i,1)=[X1{i}(freqIndexes(k));...
                X2{i}(freqIndexes(k))];
        end
    end
    Z=cell2mat(Zint');
    %periodogramme barlett
    
    for k=1:B
        C{k}=zeros(2,2);
        for i=1:4
            C{k}=C{k}+Zint{k}(2*i-1:2*i,1)*(Zint{k}(2*i-1:2*i,1))';
        end
        C{k}=C{k}/4;
    end
    
    %Calcul des crit�res
    I2=[1 0;0 1];
    c0=-2*N*B*log(pi);
    c1=c0-N*B;
    c2=c1-N*B;
    c3=c2;
    c4=c0;
    c5=c2;
    c6=c2;
    
    c0=0;
    c1=0;
    c2=0;
    c3=0;
    c4=0;
    c5=0;
    c6=0;
    
    switch(cas_de_figure)
        case 1
            for ntheta=1:Ntheta
                sum=0;
                for k=1:B
                    sum=sum+log(det(...
                        P(:,:,k,ntheta)*C{k}*(I2-P(:,:,k,ntheta))...
                        +sigma^2*(I2-P(:,:,k,ntheta)))+ ...
                        trace(  (I2-P(:,:,k,ntheta))*C{k}));
                    
                end
                J(numero_exp,ntheta)=c1-N*real(sum);
            end
        case 2
            for ntheta=1:Ntheta
                sum1=0;
                for k=1:B
                    sum1=sum1+trace(  (I2-P(:,:,k,ntheta))*C{k});
                end
                sum=0;
                for k=1:B
                    sum=sum+log(det(...
                        P(:,:,k,ntheta)*C{k}*P(:,:,k,ntheta)...
                        +(I2-P(:,:,k,ntheta))*real(sum1)/B));
                    
                end
                J(numero_exp,ntheta)=c2-N*real(sum);
            end
        case 3
            for ntheta=1:Ntheta
                sum=0;
                for k=1:B
                    sum=sum+log(det(...
                        P(:,:,k,ntheta)*C{k}*P(:,:,k,ntheta)...
                        +(I2-P(:,:,k,ntheta))*trace(  (I2-P(:,:,k,ntheta))*C{k})));
                    
                end
                J(numero_exp,ntheta)=c3-N*real(sum);
            end
        case 4
            for ntheta=1:Ntheta
%                 sum0=0;
%                 for k=1:B
%                     sum0=sum0+log(sigma^2);
%                 end
                sum1=0;
                for k=1:B
                    sum1=sum1+trace(  (I2-P(:,:,k,ntheta))*C{k})/sigma^2;
                end
%                J(numero_exp,ntheta)=c4-2*N*sum0-N*real(sum1);
                J(numero_exp,ntheta)=-N*real(sum1);
            end
        case 5
            for ntheta=1:Ntheta
                sum1=0;
                for k=1:B
                    sum1=sum1+trace(  (I2-P(:,:,k,ntheta))*C{k});
                end
                J(numero_exp,ntheta)=c5-2*N*log(real(sum1)/(2*B));
            end
        case 6
            for ntheta=1:Ntheta
                sum=0;
                for k=1:B
                    sum=sum+log(.5*trace(  (I2-P(:,:,k,ntheta))*C{k}));
                end
                J(numero_exp,ntheta)=c6-2*N*real(sum);
            end
            
    end
    
    %on prend le max
    [maxcrit_exp(numero_exp),idxmax]=max(J(numero_exp,:));
    theta_mle(numero_exp)=thetaArg(idxmax);
    fprintf('Theta mle :   %.0f  Vraie localisation:   %d\n',thetaArg(idxmax)*180/pi,round(azimuth(numero_exp)));
end
% FIN BOUCLE


%Jpos=rescale_matrix01(J);
Jpos=arrangement(J,.15);