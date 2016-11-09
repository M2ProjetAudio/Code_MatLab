clear variables
close all
clc
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
[signal_d,fs]=audioread('chasseurs.wav');
signal=signal_d(:,1)+signal_d(:,2);
%elimination des premiers echantillons nuls (le son ne commence pas
%encore); si on le fait pas => pb dans le critere qui vaudra Inf
for i=1:length(signal)
    if signal(i,1)>0
        debut=i; %enregistrement de l'indice ou le son debute
        break;
    end
end
clear signal;
signal=signal_d(debut:end,1)+signal_d(debut:end,2);

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
Taille_groupe=10000;

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

for k=freqIndexes
    for ntheta=1:Ntheta
        V=[1;H(k,ntheta)];
        P(:,:,k,ntheta)=(V*V')/(V'*V); %donnee
        %P(:,:,k,ntheta)=V/(V'*V)*V';
    end
end



% BOUCLE SUR CHAQUE POSITION for k=0 -> nb_positions
Nb_Loca=floor(length(signal)/Taille_groupe);
% azimuth a indiquer en degres
azimuth=linspace(0,360,Nb_Loca);
D=5*ones(size(azimuth));

%fenetre:
N=5+floor((Taille_groupe-taille_min_grp)/768); %hann
Taille_groupe=taille_min_grp+(N-5)*768;
% N=floor(Taille_groupe/K); %rect
B=Nbfreq;
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
for i=1:N
    X1{i}=zeros(1024,1);
    X2{i}=zeros(1024,1);
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
J6(length(D),Ntheta)=0;
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
for exp=1:length(D)
    deb=Taille_groupe*(exp-1);
    % segments extremites
    %fprintf('%d:%d\n',deb+1,deb+1024);
    x1t{1}(:,1)=x1(deb+1:deb+1024).*w;
    x2t{1}(:,1)=x2(1:1024).*w;
    X1{1}=coef*fft(x1t{1});  %zero padding
    X2{1}=coef*fft(x2t{1});
    
    for i=2:N % segments internes
        %fprintf('%d:%d\n',deb+1024*(i-1)*.75+1,deb+1024*(i-1)*.75+1024);
        x1t{i}(:,1)=x1(deb+1024*(i-1)*.75+1:deb+1024*(i-1)*.75+1024).*w;
        x2t{i}(:,1)=x2(deb+1024*(i-1)*.75+1:deb+1024*(i-1)*.75+1024).*w;
        % fft
        X1{i}=coef*fft(x1t{i});
        X2{i}=coef*fft(x2t{i});
    end
    for k=freqIndexes
        for i=1:N
            Zint{k}(2*i-1:2*i,1)=[X1{i}(k);X2{i}(k)];
        end
    end
    Z=cell2mat(Zint');
    %periodogramme barlett
    
    for k=1:B
        C{k}=zeros(2,2);
        for i=1:N
            C{k}=C{k}+Zint{k}(2*i-1:2*i,1)*(Zint{k}(2*i-1:2*i,1))';
        end
        C{k}=C{k}/N;
    end
    
    %Calcul des critères
    I2=[1 0;0 1];
    for ntheta=1:Ntheta
        sum=0;
        for k=freqIndexes
            sum=sum+log(.5*trace(  (I2-P(:,:,k,ntheta))*C{k}));
        end
        J6(exp,ntheta)=-sum;
    end
    %on prend le max
    [maxcrit_exp(exp),idxmax]=max(J6(exp,:));
    theta_mle(exp)=thetaArg(idxmax);
    disp(thetaArg(idxmax)*180/pi)
    
end
% FIN BOUCLE
mincrit=min(min(J6));
if mincrit<0
    J6pos=J6+abs(mincrit)*ones(size(J6));
    maxcrit=max(maxcrit_exp(maxcrit_exp~=Inf))+abs(mincrit);
else
    J6pos=J6;
    maxcrit=max(maxcrit_exp(maxcrit_exp~=Inf));
end

% pour eliminer les Inf (ça bousille le truc sinon)
% Les maxcrit_exp de valeur Inf venaient du fait, qu'au debut du signal
% sonore, plusieurs trames valaient 0.. (aucun son du tout)
mincrit_exp=zeros(Nb_Loca,1);
for exp=1:Nb_Loca
[maxcrit_exp(exp),idxmax]=max(J6pos(exp,:));
mincrit_exp(exp)=min(J6pos(exp,:));
J6pos(exp,:)=J6pos(exp,:)-mincrit_exp(exp);
maxcrit_exp(exp)=maxcrit_exp(exp)-mincrit_exp(exp);
end
%%
%       joue le son sur stereo
soundsc(outputSignal,fs);
%         AFFICHAGE

figure(50)
for exp=1:length(D)
    
    
    %source
    plot(D(exp)*cos(azimuth(exp)*pi/180+pi/2),...
        D(exp)*sin(azimuth(exp)*pi/180+pi/2),'sr')
    hold on
    grid on
    text(D(exp)*cos(azimuth(exp)*pi/180+pi/2),...
        D(exp)*sin(azimuth(exp)*pi/180+pi/2),'  Source sonore')
    
    
    %J6norme=(1./J6(exp,:))/(1/maxcrit)*D(exp)*.9; %norme pour la figure
    J6norme=J6pos(exp,:)/maxcrit_exp(exp)*D(exp)*.9;
    
    fill(cos(thetaArg+pi/2).*J6norme,sin(thetaArg+pi/2).*J6norme,[.8 .8 .8])
    %plot(0,0,'bo')
    circle([0,0],.5,'b','LineWidth',2)
    circle([-.5,0],.1,'b','LineWidth',2)
    circle([.5,0],.1,'b','LineWidth',2)
    text(0,-.5,'  Tete')
    line([0,0],[0,1],'LineWidth',2)
    quiver(0,0,D(exp)*cos(theta_mle(exp)+pi/2),D(exp)*sin(theta_mle(exp)+pi/2),'r')
    %quiver(0,0,0,1,'b')
    title(sprintf...
        ('version 4: Localisation critere 1:\n source en  %.0f degres et distance=%.1f m',azimuth(exp),D(exp)))
    xlabel('x')
    ylabel('y')
    axis([-D(exp)-.5 D(exp)+.5 -D(exp)-.5 D(exp)+.5])
    fig=gcf;
    fig.Color='w';
    hold off
    pause(duree_son/Nb_Loca)
end




audiowrite('chasseurs360.wav',outputSignal,fs)

