cla(handles.axes1);
cla(handles.axes2);

set(handles.axes1,'Visible','On');
set(handles.axes2,'Visible','On');

%% fig 2
axes(handles.axes2);
plot(sum(outputSignal,2));
hold on
axis off;
myStruct.barre_rouge=plot([0 0], [-1.5 1.5], 'r', 'LineWidth', 2);
%fig 1
axes(handles.axes1);
axis off;
%nbre d'experiences
nb_exp=length(D);
hold off
%carre de la source:
myStruct.carre_source=line('Visible','on','Marker','s',...
    'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','b');

hold on
%texte de la source
myStruct.texte_source=text(0.1,0.314,'source sonore','Visible','on');
%bulles des valeurs des criteres
myStruct.criteres_J_fill=fill(0,0,[.8 .8 .8],'Visible','on');
%fleche du mle
myStruct.fleche=quiver(0,0,1,1,'r','Visible','on');

%bloc de calcul
dt=length(outputSignal)/fs/nb_exp;

Jnorme{nb_exp}=zeros(size(Jpos));  myStruct.xsource(nb_exp)=0; myStruct.ysource(nb_exp)=0;
myStruct.u(nb_exp)=0; myStruct.xtext(nb_exp)=0; myStruct.ytext(nb_exp)=0;
myStruct.v(nb_exp)=0; myStruct.xfill{nb_exp}=zeros(1,360); myStruct.yfill{nb_exp}=zeros(1,360);
myStruct.xbarre_rouge(nb_exp)=0;
for num_exp=1:nb_exp
   Jnorme{num_exp}=Jpos(num_exp,:)/maxcrit_exp(num_exp)*D(num_exp)*.9;
   myStruct.xfill{num_exp}=cos(thetaArg+pi/2).*Jnorme{num_exp};
   myStruct.yfill{num_exp}=sin(thetaArg+pi/2).*Jnorme{num_exp};
   myStruct.u(num_exp)=D(num_exp)/2*cos(theta_mle(num_exp)+pi/2);
   myStruct.v(num_exp)=D(num_exp)/2*sin(theta_mle(num_exp)+pi/2);
   myStruct.xsource(num_exp)=D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2);
   myStruct.ysource(num_exp)=D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2);
   myStruct.xtext(num_exp)=(D(num_exp)+.5)*cos(azimuth(num_exp)*pi/180+pi/2);
   myStruct.ytext(num_exp)=(D(num_exp)+.5)*sin(azimuth(num_exp)*pi/180+pi/2);
   myStruct.xbarre_rouge(num_exp)=round(num_exp*dt*fs);
end

circle([0,0],.5,'b','LineWidth',2)
circle([-.5,0],.1,'b','LineWidth',2)
circle([.5,0],.1,'b','LineWidth',2)
text(0,-.5,'  Tete')
line([0,0],[0,1],'LineWidth',2)
% title(sprintf...
%     ('Resultat estimation de la Localisation.\n source en  %.0f degres et distance=%.1f m',azimuth(num_exp),D(num_exp)))
xlabel('x')
ylabel('y')
axis([-D(1)-.5 D(1)+.5 -D(1)-.5 D(1)+.5])
% fig=gcf;
% fig.Color='w';



myStruct.num_exp=1;
myStruct.nb_exp=nb_exp;

% <video>
myStruct.prendre_video=0;
if myStruct.prendre_video
   writerObj = VideoWriter('vid1','MPEG-4');
   writerObj.FrameRate = evalin('base','fps');
   myStruct.frame{myStruct.nb_exp}=getframe(Menu2);
   for k=1:myStruct.nb_exp
       myStruct.frame{k}=myStruct.frame{myStruct.nb_exp};
   end
    open(writerObj);
end
% </video>


p=audioplayer(outputSignal,fs);
set(p, 'UserData', myStruct);
set(p, 'TimerFcn', @synchroCallback);
set(p,'StartFcn',@synchroPremier);
set(p, 'TimerPeriod',dt);


playblocking(p)

myStruct=get(p, 'UserData');
if myStruct.prendre_video
    for k=1:nb_exp
         writeVideo(writerObj,myStruct.frame{k});
    end
    
    close(writerObj);
end



if parole
    audiowrite([wav_name,'360.wav'],outputSignal,fs)
elseif gauss
    audiowrite('gauss360.wav',outputSignal,fs)
end
