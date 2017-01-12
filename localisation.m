function sortie_loca=localisation(gauss,signal,wav_name,cas_de_figure,duree_son)

fps=evalin('base','fps');
fs=44100;

sortie_space=spatialisation(gauss,signal,wav_name,cas_de_figure,duree_son);
Lframe=sortie_space.Lframe;
Ntheta=sortie_space.Ntheta;
B=sortie_space.B;
P=sortie_space.P;
thetaArg=sortie_space.thetaArg;
freqIndexes=sortie_space.freqIndexes;



% BOUCLE SUR CHAQUE POSITION for k=0 -> nb_positions
%Nb_Loca=floor(length(signal)/Taille_de_1_position);
Nb_Loca=floor(duree_son*fps);
while (length(signal)/Nb_Loca)< Lframe*2.5
    fps=fps-1;
    Nb_Loca=floor(duree_son*fps);
    fprintf('!!!!!!!!!!! DIMINUTION FPS !!!!!!!!!!\n\n');
end
fprintf('!!!!!!!!!!! FPS  vaut  %d !!!!!!!!!!\n\n',fps);
assignin('base','fps',fps);
fprintf('Il y aura %d positionnements et autant de localisations\n',Nb_Loca);
% azimuth a indiquer en degres
azimuth=linspace(0,360,Nb_Loca);
D=5*ones(size(azimuth));

Taille_de_1_position=floor(length(signal)/length(D));
%fenetre:
N=4; %5 hann par localisation
%Taille_de_1_position=20000;
% N=floor(Taille_de_1_position/K); %rect
%B=Nbfreq;
%Allocations de memoire
J1(Ntheta,1)=0;
x1t{N}=NaN;
x2t{N}=NaN;
for i=1:N
    x1t{i}=zeros(Lframe,1);
    x2t{i}=zeros(Lframe,1);
end
X1{N}=NaN;
X2{N}=NaN;
for i=1:4
    X1{i}=zeros(Lframe,1);
    X2{i}=zeros(Lframe,1);
end
C{B}=NaN;
for k=1:B
    C{k}=zeros(2,2);
end
Zint{B}=NaN;
for k=1:B
    Zint{k}=zeros(2*N,1);
end

outputSignal=zeros(Taille_de_1_position*length(D),2);
J(length(D),Ntheta)=0;
theta_mle=zeros(size(D));
maxcrit_exp=zeros(size(D));
disp('Localisation estimee:\n')

hrir = simulator.DirectionalIR(xml.dbGetFile...
    ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa'));

for numero_exp=1:length(D)
    
    impulseResponse = hrir.getImpulseResponses(azimuth(numero_exp));
    %fprintf('%d:%d\n',(numero_exp-1)*Taille_de_1_position+1,numero_exp*Taille_de_1_position)
    signal_tr=signal((numero_exp-1)*Taille_de_1_position+1:numero_exp*Taille_de_1_position);
    % application de la hrir aux 2 oreilles
    left_ear=conv(signal_tr,impulseResponse.left);
    right_ear=conv(signal_tr,impulseResponse.right);
    
    outputSignal((numero_exp-1)*Taille_de_1_position+1:...
        numero_exp*Taille_de_1_position,:) = ...
        [left_ear(1:Taille_de_1_position) ...
        right_ear(1:Taille_de_1_position)];
end
%%
% Ajout de bruit sur le signal d'entree
Psig=mean(mean(outputSignal,2).^2);
% if Psig<1
%     sigma=Psig/10;
% else
    sigma=sqrt(Psig/10)/10
% end

bruit1=randn(size(outputSignal,1),1);
bruit2=randn(size(outputSignal,1),1);

[b,a]=butter(3,5000/(fs/2),'low');
bruit1=filter(b,a,bruit1);
bruit2=filter(b,a,bruit2);

outputSignal=outputSignal+sigma*[bruit1,bruit2];
%re-scale du outputSignal entre -1 et 1
outputSignal=rechelonner(outputSignal);


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


Jpos=arrangement(J,.15);


  sortie_loca.outputSignal=outputSignal;
  sortie_loca.D=D;
  sortie_loca.sigma=sigma;
  sortie_loca.azimuth=azimuth;
  sortie_loca.Jpos=Jpos;
  sortie_loca.maxcrit_exp=maxcrit_exp;
  sortie_loca.theta_mle=theta_mle;
  sortie_loca.Nb_Loca=Nb_Loca;
  
  




  
 
