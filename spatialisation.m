function sortie_space=spatialisation(gauss,signal,wav_name,cas_de_figure,duree_son)
fps=evalin('base','fps');
fs=44100;



% GET HRTF from HRIR AND INIT LOCALISATION

hrtffile = xml.dbGetFile(...
    'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa');
hrir_SOFA= SOFAload(hrtffile);

Lframe = 512; % taille d'une trame
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






  sortie_space.Lframe=Lframe;
  sortie_space.B=B;
  sortie_space.Ntheta=Ntheta;
  sortie_space.P=P;
  sortie_space.thetaArg=thetaArg;
  sortie_space.freqIndexes=freqIndexes;

  
  
  


