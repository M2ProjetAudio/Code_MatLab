
% var a importer du menu%
%cas_de_figure=6;
%wav_name='chasseurs';
%  debut algo %

% ATTENTION : ALGORITHME TRES LENT
% CONSEIL: LANCER le premier bloc (ctrl+entree) pour execution des
% localisations
% Pui lancer second bloc (ctrl+entree) apres les '%%' pour ecoute et
% affichage simultanee du resultat

% ATTENTION: avoir le fichier circle.m dans le mm repertoire pour le dessin


%   test spatialisation d'un signal qui fait un 360 autour de l'auditeur.
% un signal sonore de 25 sec
% s'ecoute en stereo.
% je voulais:
% Resultat: fonctionne comme attendu !

% INIT SOURCE -> c'est le fichier son non-spatialise
% fs=44100;
% signal=randn(fs,1);

if parole
    if strcmp(wav_name(end-4:end),'.wav')
        file=wav_name;
        wav_name(end-4:end)=[];
    else
        file=[wav_name, '.wav'];
    end
    
    
    [signal_d,fs]=audioread(file);
    if size(signal_d,2)==2
        signal=signal_d(:,1)+signal_d(:,2);
    else
        signal=signal_d;
    end
    %elimination des premiers echantillons nuls (le son ne commence pas
    %encore); si on le fait pas => pb dans le critere qui vaudra Inf
    for i=1:length(signal)
        if signal(i,1)>0
            debut=i; %enregistrement de l'indice ou le son debute
            break;
        end
    end
    clear signal;
    if size(signal_d,2)==2
    signal=signal_d(debut:end,1)+signal_d(debut:end,2);
    else signal=signal_d(debut:end,1);
    end
    
elseif gauss
    fs=44100;
    signal=randn(fs*8,1);
end
duree_son = length(signal) / fs;% => 25 sec

% o=audioplayer(signal,fs);
% pour jouer le fichier    >>  play(o)

% GET HRTF from HRIR AND INIT LOCALISATION

hrtffile = xml.dbGetFile(...
    'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa');
hrir_SOFA= SOFAload(hrtffile);

Lframe = 1024; % ??
az=(-180:179);


HRIR = hrir_SOFA.Data.IR;
HRIR = permute(HRIR,[3 2 1]);
Nangles = length(az);


%taille de chaque groupe (nb de points )
taille_min_grp=4096; %=> correspond à 5trames de 1024 ech
%Taille_groupe=4096+2*768;


K=1024; % nombre de points par trame
% Rq: il y aura taille_groupe/K trames dans chaque groupe de localisation
% si la fenetre est rect
% Plus de trames si fenetres de Hann (a cause du recouvrement)
Noverlap=floor(K/2); % recouvrement des trames
%nombre de trames
Nf=4;
Fmax = 40000;
Nfreqs = fix(Fmax*Lframe/fs) + 1;

maxIndex=floor(K/2)+1;
B=maxIndex;
freqIndexes=1:maxIndex;
Nbfreq=length(freqIndexes);
% grille des azimuths testes
thetaArg=az*pi/180;
thetaArg2 = [-185*pi/180 -180*pi/180 thetaArg 185*pi/180];
% nombre d'azimuths testes
Ntheta=length(thetaArg);


HRTF = NaN(Lframe,2,360);
% compute the anechoic hrtfs by fourier transforming irs. zero padd irs if
% irs length is smaller than a frame length, otherwise truncate irs.

if Lframe <= 2048
    for i=1:length(az)
        HRTF(:,:,i)=fft(HRIR(1:Lframe,:,i));
    end
else
    for i=1:length(az)
        HRTF(:,:,i)=fft([HRIR(:,:,i) ; zeros(Lframe-2048,2)]);
    end
end

H = squeeze(HRTF(:,2,:)./HRTF(:,1,:));

% calcul du projecteur pour la localisation
P=NaN(2,2,Nbfreq,Ntheta);

for k=1:B
    for ntheta=1:Ntheta
        V=[1;H(k,ntheta)];
        P(:,:,k,ntheta)=(V*V')/(V'*V); %donnee
        %P(:,:,k,ntheta)=V/(V'*V)*V';
    end
end


% BOUCLE SUR CHAQUE POSITION for k=0 -> nb_positions
%Nb_Loca=floor(length(signal)/Taille_groupe);
Nb_Loca=72;
% azimuth a indiquer en degres
azimuth=linspace(0,360,Nb_Loca);
D=5*ones(size(azimuth));

Taille_groupe=floor(length(signal)/length(D));
%fenetre:
N=4; %5 hann par localisation
%Taille_groupe=20000;
% N=floor(Taille_groupe/K); %rect
B=Nbfreq;
B=128;
%Allocations de memoire
J1(Ntheta,1)=0;
x1t{N}=NaN;
x2t{N}=NaN;
for i=1:N
    x1t{i}=zeros(1024,1);
    x2t{i}=zeros(1024,1);
end
X1{N}=NaN;
X2{N}=NaN;
for i=1:4
    X1{i}=zeros(128,1);
    X2{i}=zeros(128,1);
end
C{B}=NaN;
for k=1:B
    C{k}=zeros(2,2);
end
Zint{B}=NaN;
for k=1:B
    Zint{k}=zeros(2*N,1);
end

outputSignal=zeros(Taille_groupe*length(D),2);
J(length(D),Ntheta)=0;
theta_mle=zeros(size(D));
maxcrit_exp=zeros(size(D));
disp('Localisation estimee:\n')

hrir = simulator.DirectionalIR(xml.dbGetFile...
    ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa'));

for exp=1:length(D)
    % SI k=0
    %AFFICHAGE DE LA SITUATION INITIALE
    ...
        %  SINON
    %CALCUL DE LA COMMANDE -> u_k = u_0 pour la BO
    ...
        %MAJ POSITION
    ...
        %PREPARATION DE LA NOUVELLE MESURE -> SPATIALISATION DU SIGNAL
    % D(exp) -> distance source - capteur
    %
    %     if D(exp) >= 3
    %         hrir = simulator.DirectionalIR(xml.dbGetFile...
    %             ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa'));
    %     elseif (D(exp) <3 && D(exp)>=2)
    %         hrir = simulator.DirectionalIR(xml.dbGetFile...
    %             ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_2m.sofa'));
    %     else
    %         hrir = simulator.DirectionalIR(xml.dbGetFile...
    %             ('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_1m.sofa'));
    %     end
    %
    
    impulseResponse = hrir.getImpulseResponses(azimuth(exp));
    
    %fprintf('%d:%d\n',(exp-1)*Taille_groupe+1,exp*Taille_groupe)
    signal_tr=signal((exp-1)*Taille_groupe+1:exp*Taille_groupe);
    outputSignal((exp-1)*Taille_groupe+1:exp*Taille_groupe,:) = [conv(signal_tr,impulseResponse.left,'same') ...
        conv(signal_tr,impulseResponse.right,'same')];
end
%%
% Ajout de bruit sur le signal d'entree
Psig=mean(mean(outputSignal.^2));
if Psig<1
    sigma=Psig/10;
else
    sigma=sqrt(Psig/10);
end
outputSignal=outputSignal+sigma*randn(size(outputSignal));
x1=outputSignal(:,1);
x2=outputSignal(:,2);
% decoupage en N intervalles de taille 1024
w=hanning(1024);
%w=1;
coef=sqrt(N/duree_son);
%%

for exp=1:length(D)-2
    deb=Taille_groupe*(exp-1)+15+Taille_groupe/2; %debut a peu pres au milieu dun sig localise
    % segments extremites
    fprintf('deb: %d\n',deb)
    %fprintf('%d:%d\n',deb+1,deb+1024);
    x1t{1}(:,1)=x1(deb+1:deb+1024).*w;
    x2t{1}(:,1)=x2(deb+1:deb+1024).*w;
    X1{1}=coef*fft(x1t{1});  
    X2{1}=coef*fft(x2t{1}); 
    
    for i=2:4 % segments internes
        %fprintf('%d:%d\n',deb+1024*(i-1)*.75+1,deb+1024*(i-1)*.75+1024);
        x1t{i}(:,1)=x1(deb+512*(i-1)+1:deb+512*(i-1)+1024).*w;
        x2t{i}(:,1)=x2(deb+512*(i-1)+1:deb+512*(i-1)+1024).*w;
        % fft
        X1{i}=coef*fft(x1t{i});
        X2{i}=coef*fft(x2t{i});
    end
    for k=1:B
        for i=1:4
            Zint{k}(2*i-1:2*i,1)=[X1{i}(k);X2{i}(k)];
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
    
    %Calcul des critères
    I2=[1 0;0 1];
    c0=-2*5*B*log(pi);
    c1=c0-5*B;
    c2=c1-5*B;
    c3=c2;
    c4=c0;
    c5=c2;
    c6=c2;
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
                J(exp,ntheta)=c1-5*real(sum);
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
                J(exp,ntheta)=c2-5*real(sum);
            end
        case 3
            for ntheta=1:Ntheta
                sum=0;
                for k=1:B
                    sum=sum+log(det(...
                        P(:,:,k,ntheta)*C{k}*P(:,:,k,ntheta)...
                        +(I2-P(:,:,k,ntheta))*trace(  (I2-P(:,:,k,ntheta))*C{k})));
                    
                end
                J(exp,ntheta)=c3-5*real(sum);
            end
        case 4
            for ntheta=1:Ntheta
                sum0=0;
                for k=1:B
                    sum0=sum0+log(sigma^2);
                end
                sum1=0;
                for k=1:B
                    sum1=sum1+trace(  (I2-P(:,:,k,ntheta))*C{k})/sigma^2;
                end
                J(exp,ntheta)=c4-2*5*sum0-5*real(sum1);
            end
        case 5
            for ntheta=1:Ntheta
                sum1=0;
                for k=1:B
                    sum1=sum1+trace(  (I2-P(:,:,k,ntheta))*C{k});
                end
                J(exp,ntheta)=c5-2*5*log(real(sum1)/(2*B));
            end
        case 6
            for ntheta=1:Ntheta
                sum=0;
                for k=1:B
                    sum=sum+log(.5*trace(  (I2-P(:,:,k,ntheta))*C{k}));
                end
                J(exp,ntheta)=c6-2*5*real(sum);
            end
            
    end
    
    %on prend le max
    [maxcrit_exp(exp),idxmax]=max(J(exp,:));
    theta_mle(exp)=thetaArg(idxmax);
    fprintf('Theta mle :')
    disp(thetaArg(idxmax)*180/pi)
     fprintf('Vraie localisation:  %d\n',round(azimuth(exp)));
    
end
% FIN BOUCLE
mincrit=min(min(J));
if mincrit<0
    Jpos=J+abs(mincrit)*ones(size(J));
    maxcrit=max(maxcrit_exp(maxcrit_exp~=Inf))+abs(mincrit);
else
    Jpos=J;
    maxcrit=max(maxcrit_exp(maxcrit_exp~=Inf));
end

% pour eliminer les Inf (ça bousille le truc sinon)
% Les maxcrit_exp de valeur Inf venaient du fait, qu'au debut du signal
% sonore, plusieurs trames valaient 0.. (aucun son du tout)
mincrit_exp=zeros(Nb_Loca,1);
for exp=1:Nb_Loca
    [maxcrit_exp(exp),idxmax]=max(Jpos(exp,:));
    mincrit_exp(exp)=min(Jpos(exp,:));
    Jpos(exp,:)=Jpos(exp,:)-mincrit_exp(exp);
    maxcrit_exp(exp)=maxcrit_exp(exp)-mincrit_exp(exp);
end