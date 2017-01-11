function sortie_space=spatialisation(gauss,signal,wav_name)
fps=evalin('base','fps');
fs=44100;
duree_son = length(signal) / fs;% => 25 sec

% GET HRTF from HRIR AND INIT LOCALISATION

hrtffile = xml.dbGetFile(...
    'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa');
hrir_SOFA= SOFAload(hrtffile);

Lframe = 128; % taille d'une trame
az=(-180:179);


HRIR = hrir_SOFA.Data.IR;
HRIR = permute(HRIR,[3 2 1]);
Nangles = length(az);


 % nombre de points par trame
% Rq: il y aura taille_groupe/K trames dans chaque groupe de localisation
% si la fenetre est rect
% Plus de trames si fenetres de Hann (a cause du recouvrement)
 % recouvrement des trames
%nombre de trames
Nf=4;
Fmax = 40000;
Nfreqs = fix(Fmax*Lframe/fs) + 1;

%maxIndex=floor(K/2)+1;
B=128;
freqIndexes=round(linspace(1,round(Lframe) ,B));
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
        V=[1;H(freqIndexes(k),ntheta)];
        P(:,:,k,ntheta)=(V*V')/(V'*V); %donnee
        %P(:,:,k,ntheta)=V/(V'*V)*V';
    end
end


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





  sortie_space.N=N;
  sortie_space.outputSignal=outputSignal;
  sortie_space.D=D;
  sortie_space.Lframe=Lframe;
  sortie_space.B=B;
  sortie_space.Ntheta=Ntheta;
  sortie_space.P=P;
  sortie_space.sigma=sigma;
  sortie_space.thetaArg=thetaArg;
  sortie_space.azimuth=azimuth;
  sortie_space.Taille_de_1_position=Taille_de_1_position;
  sortie_space.freqIndexes=freqIndexes;
  