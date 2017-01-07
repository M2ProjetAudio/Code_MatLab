cla(handles.axes1);
%nbre d'experiences
nb_exp=length(D);
hold off
%carre de la source:
myStruct.carre_source=line('Visible','on');

hold on
%texte de la source
myStruct.texte_source=text(0.1,0.314,'source sonore','Visible','on');
%bulles des valeurs des criteres
myStruct.criteres_J_fill=fill(0,0,[.8 .8 .8],'Visible','on');
%fleche du mle
myStruct.fleche=quiver(0,0,1,1,'r','Visible','on');

Jnorme{nb_exp}=zeros(size(Jpos));  myStruct.xsource(nb_exp)=0; myStruct.ysource(nb_exp)=0;
myStruct.u(nb_exp)=0; myStruct.xtext(nb_exp)=0; myStruct.ytext(nb_exp)=0;
myStruct.v(nb_exp)=0; myStruct.xfill{nb_exp}=zeros(1,360); myStruct.yfill{nb_exp}=zeros(1,360);
for num_exp=1:nb_exp
   Jnorme{num_exp}=Jpos(num_exp,:)/maxcrit_exp(num_exp)*D(num_exp)*.9;
   myStruct.xfill{num_exp}=cos(thetaArg+pi/2).*Jnorme{num_exp};
   myStruct.yfill{num_exp}=sin(thetaArg+pi/2).*Jnorme{num_exp};
   myStruct.u(num_exp)=D(num_exp)/2*cos(theta_mle(num_exp)+pi/2);
   myStruct.v(num_exp)=D(num_exp)/2*sin(theta_mle(num_exp)+pi/2);
   myStruct.xsource(num_exp)=D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2);
   myStruct.ysource(num_exp)=D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2);
   myStruct.xtext(num_exp)=D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2);
   myStruct.ytext(num_exp)=D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2);
end

circle([0,0],.5,'b','LineWidth',2)
circle([-.5,0],.1,'b','LineWidth',2)
circle([.5,0],.1,'b','LineWidth',2)
text(0,-.5,'  Tete')
line([0,0],[0,1],'LineWidth',2)
title(sprintf...
    ('Resultat estimation de la Localisation.\n source en  %.0f degres et distance=%.1f m',azimuth(num_exp),D(num_exp)))
xlabel('x')
ylabel('y')
axis([-D(1)-.5 D(1)+.5 -D(1)-.5 D(1)+.5])
fig=gcf;
fig.Color='w';


dt=length(outputSignal)/fs/nb_exp;
myStruct.num_exp=1;

p=audioplayer(outputSignal,fs);
set(p, 'UserData', myStruct);
set(p, 'TimerFcn', @synchroCallback);
set(p, 'TimerPeriod',dt);
playblocking(p)

% %         AFFICHAGE
% for num_exp=1:nb_exp
%    
%     carre_source.XData=xsource(num_exp);
%     carre_source.YData=ysource(num_exp);
%     texte_source.Position=[xtext(num_exp) ytext(num_exp) 0];
%     fleche.UData=u(num_exp);
%     fleche.VData=v(num_exp);
%     criteres_J_fill.XData=xfill{num_exp};
%     criteres_J_fill.YData=yfill{num_exp};
%     
%     if num_exp==1
%         carre_source=line('Visible','on');
%         %texte de la source
%         texte_source=text(0,0,'source sonore','Visible','on');
%         %bulles des valeurs des criteres
%         criteres_J_fill=fill(0,0,[.8 .8 .8],'Visible','on');
%         %fleche du mle
%         fleche=quiver(0,0,1,1,'r','Visible','on');
%     end
%     
%     pause(dt)
% end
% 
% 


if parole
    audiowrite([wav_name,'360.wav'],outputSignal,fs)
elseif gauss
    audiowrite('gauss360.wav',outputSignal,fs)
end
