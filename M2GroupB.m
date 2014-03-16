
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Group 'B'  Class: M2                   %
%                                                   %
%              Mickdhad Bava 11M202                 %
%              Dhruv Chand   11M204                 %
%              Safwan CH     11M237                 %
%              Vinit Vardhan 11M271                 %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;
clc;
close all;

[Selection,ok] = listdlg('ListString',{'Graphical Analysis' 'Static Force Analysis' 'Dynamic Force Analysis' 'Alternative Simulation','Inverse Static Force Analysis','Testing Mode'},'SelectionMode','single','Name','What do you want to do?','ListSize',[300 150]);
if(ok==0)
    msgbox('Bye!');

else
      if(Selection==6  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
display('Analysis of a Four Bar Mechanism')
display('All angular variables are positive in anticlockwise sense.')
display('Default value of 10 will be taken in the case of no input.')

input = inputdlg({'Enter length of link 1 (m)','Enter length of link 2 (m)','Enter length of link 3 (m)','Enter length of link 4 (m)','Enter the rpm of the crank','Enter the angular acceleration of the crank(rad/s^2)','Enter theta value'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1','1'});
a= str2num(input{2});
if isempty(a)
    a = 50;
end


b=str2num(input{3});
if isempty(b)
    b = 60;
end

c=str2num(input{4});
if isempty(c)
    c = 50;
end

d=str2num(input{1});
if isempty(d)
    d = 50;
end

w2=str2num(input{5});
if isempty(w2)
    w2 = 1;
end

e2=str2num(input{6});
if isempty(e2)
    e2 = 10;
end

t2=deg2rad(str2num(input{7}));
if isempty(t2)
    t2 = 10;
end

w2=2*pi*w2/60;

if(d>(a+b+c))
   msgbox('Incompatible link lengths!')    
else

    
    %Position Analysis
    z2 =a*exp(1i*t2);
    k1=d/a;
    k2=d/c;
    k3=(a^2+c^2+d^2-b^2)/(2*a*c);
    k4=d/b;
    k5=((c^2-d^2-a^2-b^2)/(2*a*b));
    A=cos(t2)-k1-(k2*cos(t2))+k3;
    B=(-2*sin(t2));
    C=k1-((k2+1)*cos(t2))+k3;
    D=cos(t2)-k1+(k4*cos(t2))+k5;
    E=(-2*sin(t2));
    F=k1+((k4-1)*cos(t2))+k5;
    if((B^2-(4*A*C))<0)
       
    else
        if((a==b&&b==c&&c==d)||(a==c&&b==d))
            t4=t2;
        else
    t4=2*atan((-B-sqrt(B^2-(4*A*C)))/(2*A));%Computing angle theta4
        end
    t3=2*atan((-E-sqrt(E^2-(4*D*F)))/(2*D));%Computing angle theta3
    z4=c*exp(1i*t4);
   
    %Velocity Computation
    w3=(a*w2*sin(t4-t2))/(b*sin(t3-t4));%angular velocity of link AB
    w4=(a*w2*sin(t2-t3))/(c*sin(t4-t3));% angular velocity of link O4B
   r2=a*w2;%magnitudes of velocity vector V2
   r3=b*w3;%magnitudes of velocity vector V3
   r4=c*w4;%magnitudes of velocity vector V4
    V2=r2*exp(1i*((pi/2)+t2));%Representing in complex form 
    V3=r3*exp(1i*((pi/2)+t3));%Representing in complex form 
    V4=r4*exp(1i*((pi/2)+t4));%Representing in complex form     
    %Acceleration Computation 
    P=c*sin(t4);
    Q=b*sin(t3);
    R=a*e2*sin(t2)+(a*w2^2*cos(t2))+(b*w3^2*cos(t3))-(c*w4^2*cos(t4));
    S=c*cos(t4);
    T=b*cos(t3);
    U=a*e2*cos(t2)-(a*w2^2*sin(t2))-(b*w3^2*sin(t3))+(c*w4^2*sin(t4));
    e3=((R*S-P*U)/(P*T-Q*S));%alpha 3
    e4=((R*T-Q*U)/(P*T-Q*S));%alpha 4
  
    %Acceleration A2
    A2t=a*e2*exp(1i*(t2+(pi/2)));%tangential component
    A2n=a*(w2)^2*exp(1i*(t2+pi));%radial component
    A2=A2t+A2n;
    %Acceleration A3
    A3t=b*e3*exp(1i*((pi/2)+t3)) ;%tangential component 
    A3n=b*(w3)^2*exp(1i*(t3+pi)) ;%radial component
    A3=A3t+A3n;
    %Acceleration A4
    A4t=c*e4*exp(1i*((pi/2)+t4)) ;%tangential component
    A4n=c*(w4)^2*exp(1i*(t4+pi)) ;%radial component
     A4=A4t+A4n;
     
     %%transmission angle
     t_mu=acos((a^2+d^2-b^2-c^2-2*a*d*cos(t2))/(-2*b*c))
     %%%%%%%%Plotting Code%%%%%%%%%
     
    %Plotting Position Vector Configuration Diagram
    subplot(2,2,1), plot([0 real(z2)],[0 imag(z2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
    hold on;
    axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    plot([real(z2) (d+real(z4))],[imag(z2) imag(z4)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold on;
    plot([(d+real(z4)) d],[imag(z4) 0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold off;
    title('Simulation') ;
    text(0,-0.1,'O2','fontsize',10);%naming joints
    text(real(z2),imag(z2),'A','fontsize',10);
    text(d+real(z4),imag(z4),'B','fontsize',10);
    text(d,-0.1,'O4','fontsize',10);
    
      %Plotting Acceleration Vectors
    subplot(2,2,4),quiver(0,0,real(A2n),imag(A2n));    hold on;
    quiver(real(A2n),imag(A2n),1.1*(real(A2t)),(1.1*(imag(A2t))));    hold on ;
    quiver(0,0,1.1*(real(A2)),1.1*(imag(A2))) ;    hold on ;
    quiver(0,0,real(A4n),imag(A4n));    hold on;
    quiver(real(A4n),imag(A4n),1.1*(real(A4t)),1.1*(imag(A4t)));    hold on;
    quiver(0,0,real(A4),imag(A4));    hold on;
    quiver(real(A2),imag(A2),real(A3n),imag(A3n));    hold on;
    quiver(real(A2),imag(A2),real(A3),imag(A3));    hold on;
    quiver((real(A2)+real(A3n)),(imag(A2)+imag(A3n)),real(A3t),imag(A3t));    hold on;
    axis manual;   axis([-a-1 d+c+1 -(c+a) c+a])  ; hold off;
    title('Acceleration Vector Diagram:') ;
    
    text(0,0,'O','fontsize',10);
    text(real(A2n)+0.1,1.1*imag(A2n)+0.1,'A1','fontsize',8);
    text(real(A2)+0.1,1.1*imag(A2)+0.1,'A','fontsize',8);
    text(real(A4n)-0.1,1.1*imag(A4n)-0.1,'B1','fontsize',8);
    text(real(A4)-0.1,1.1*imag(A4)-0.1,'B','fontsize',8);
   
    
    %Plotting Velocity Vectors
    subplot(2,2,3),quiver(0,0,real(V2),imag(V2),0);    hold on;
      axis manual;   axis([(-a-1)/5 (d+c+1)/5 -(c+a)/5 (c+a)/6])  ;
    quiver(0,0,real(V4),imag(V4),0);
    hold on;     
    quiver(real(V2),imag(V2),(real(V4)-real(V2)),(imag(V4)-imag(V2)),0);    hold on;
    hold off;
    title('Velocity Vector Diagram') ;
    text(0,0,'o','fontsize',10);
    text(real(V2)+0.01,imag(V2)+0.01*imag(V2),'a','fontsize',8);
    text(real(V4)-0.01,imag(V4)-0.01*imag(V4),'b','fontsize',8); 
  
    subplot(2,2,2),
    axis manual;
    axis([0,2*pi,0,pi]);
    plot(t2,t_mu,'-rs','LineWidth',1,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',2);
     title('Transmission Angle v/s Input Angle Plot.') ;
   hold on;
    pause(5)
    
    display(d+real(z4));
     display(imag(z4));
    %dipsplay all necessary vars
    
    
    
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
    
      
    if(Selection==1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
display('Analysis of a Four Bar Mechanism')
display('All angular variables are positive in anticlockwise sense.')
display('Default value of 10 will be taken in the case of no input.')

input = inputdlg({'Enter length of link 1 (m)','Enter length of link 2 (m)','Enter length of link 3 (m)','Enter length of link 4 (m)','Enter the rpm of the crank','Enter the angular acceleration of the crank(rad/s^2)'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});
a= str2num(input{2});
if isempty(a)
    a = 10;
end


b=str2num(input{3});
if isempty(b)
    b = 10;
end

c=str2num(input{4});
if isempty(c)
    c = 10;
end

d=str2num(input{1});
if isempty(d)
    d = 10;
end

w2=str2num(input{5});
if isempty(w2)
    w2 = 10;
end

e2=str2num(input{6});
if isempty(e2)
    e2 = 10;
end

w2=2*pi*w2/60;

if(d>(a+b+c))
   msgbox('Incompatible link lengths!')    
else
for t2=0:pi/50:2*pi
    
    %Position Analysis
    z2 =a*exp(1i*t2);
    k1=d/a;
    k2=d/c;
    k3=(a^2+c^2+d^2-b^2)/(2*a*c);
    k4=d/b;
    k5=((c^2-d^2-a^2-b^2)/(2*a*b));
    A=cos(t2)-k1-(k2*cos(t2))+k3;
    B=(-2*sin(t2));
    C=k1-((k2+1)*cos(t2))+k3;
    D=cos(t2)-k1+(k4*cos(t2))+k5;
    E=(-2*sin(t2));
    F=k1+((k4-1)*cos(t2))+k5;
    if((B^2-(4*A*C))<0)
       
    else
        if((a==b&&b==c&&c==d)||(a==c&&b==d))
            t4=t2;
        else
    t4=2*atan((-B-sqrt(B^2-(4*A*C)))/(2*A));%Computing angle theta4
        end
    t3=2*atan((-E-sqrt(E^2-(4*D*F)))/(2*D));%Computing angle theta3
    z4=c*exp(1i*t4);
   
    %Velocity Computation
    w3=(a*w2*sin(t4-t2))/(b*sin(t3-t4));%angular velocity of link AB
    w4=(a*w2*sin(t2-t3))/(c*sin(t4-t3));% angular velocity of link O4B
   r2=a*w2;%magnitudes of velocity vector V2
   r3=b*w3;%magnitudes of velocity vector V3
   r4=c*w4;%magnitudes of velocity vector V4
    V2=r2*exp(1i*((pi/2)+t2));%Representing in complex form 
    V3=r3*exp(1i*((pi/2)+t3));%Representing in complex form 
    V4=r4*exp(1i*((pi/2)+t4));%Representing in complex form     
    %Acceleration Computation 
    P=c*sin(t4);
    Q=b*sin(t3);
    R=a*e2*sin(t2)+(a*w2^2*cos(t2))+(b*w3^2*cos(t3))-(c*w4^2*cos(t4));
    S=c*cos(t4);
    T=b*cos(t3);
    U=a*e2*cos(t2)-(a*w2^2*sin(t2))-(b*w3^2*sin(t3))+(c*w4^2*sin(t4));
    e3=((R*S-P*U)/(P*T-Q*S));%alpha 3
    e4=((R*T-Q*U)/(P*T-Q*S));%alpha 4
  
    %Acceleration A2
    A2t=a*e2*exp(1i*(t2+(pi/2)));%tangential component
    A2n=a*(w2)^2*exp(1i*(t2+pi));%radial component
    A2=A2t+A2n;
    %Acceleration A3
    A3t=b*e3*exp(1i*((pi/2)+t3)) ;%tangential component 
    A3n=b*(w3)^2*exp(1i*(t3+pi)) ;%radial component
    A3=A3t+A3n;
    %Acceleration A4
    A4t=c*e4*exp(1i*((pi/2)+t4)) ;%tangential component
    A4n=c*(w4)^2*exp(1i*(t4+pi)) ;%radial component
     A4=A4t+A4n;
     
     %%transmission angle
     t_mu=acos((a^2+d^2-b^2-c^2-2*a*d*cos(t2))/(-2*b*c))
     %%%%%%%%Plotting Code%%%%%%%%%
     
    %Plotting Position Vector Configuration Diagram
    subplot(2,2,1), plot([0 real(z2)],[0 imag(z2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
    hold on;
    axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    plot([real(z2) (d+real(z4))],[imag(z2) imag(z4)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold on;
    plot([(d+real(z4)) d],[imag(z4) 0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold off;
    title('Simulation') ;
    text(0,-0.1,'O2','fontsize',10);%naming joints
    text(real(z2),imag(z2),'A','fontsize',10);
    text(d+real(z4),imag(z4),'B','fontsize',10);
    text(d,-0.1,'O4','fontsize',10);
    
      %Plotting Acceleration Vectors
    subplot(2,2,4),quiver(0,0,real(A2n),imag(A2n));    hold on;
    quiver(real(A2n),imag(A2n),1.1*(real(A2t)),(1.1*(imag(A2t))));    hold on ;
    quiver(0,0,1.1*(real(A2)),1.1*(imag(A2))) ;    hold on ;
    quiver(0,0,real(A4n),imag(A4n));    hold on;
    quiver(real(A4n),imag(A4n),1.1*(real(A4t)),1.1*(imag(A4t)));    hold on;
    quiver(0,0,real(A4),imag(A4));    hold on;
    quiver(real(A2),imag(A2),real(A3n),imag(A3n));    hold on;
    quiver(real(A2),imag(A2),real(A3),imag(A3));    hold on;
    quiver((real(A2)+real(A3n)),(imag(A2)+imag(A3n)),real(A3t),imag(A3t));    hold on;
    axis manual;   axis([-a-1 d+c+1 -(c+a) c+a])  ; hold off;
    title('Acceleration Vector Diagram:') ;
    
    text(0,0,'O','fontsize',10);
    text(real(A2n)+0.1,1.1*imag(A2n)+0.1,'A1','fontsize',8);
    text(real(A2)+0.1,1.1*imag(A2)+0.1,'A','fontsize',8);
    text(real(A4n)-0.1,1.1*imag(A4n)-0.1,'B1','fontsize',8);
    text(real(A4)-0.1,1.1*imag(A4)-0.1,'B','fontsize',8);
   
    
    %Plotting Velocity Vectors
    subplot(2,2,3),quiver(0,0,real(V2),imag(V2),0);    hold on;
      axis manual;   axis([(-a-1)/5 (d+c+1)/5 -(c+a)/5 (c+a)/6])  ;
    quiver(0,0,real(V4),imag(V4),0);
    hold on;     
    quiver(real(V2),imag(V2),(real(V4)-real(V2)),(imag(V4)-imag(V2)),0);    hold on;
    hold off;
    title('Velocity Vector Diagram') ;
    text(0,0,'o','fontsize',10);
    text(real(V2)+0.01,imag(V2)+0.01*imag(V2),'a','fontsize',8);
    text(real(V4)-0.01,imag(V4)-0.01*imag(V4),'b','fontsize',8); 
  
    subplot(2,2,2),
    axis manual;
    axis([0,2*pi,0,pi]);
    plot(t2,t_mu,'-rs','LineWidth',1,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',2);
     title('Transmission Angle v/s Input Angle Plot.') ;
   hold on;
    pause(0.001)
    end
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
if(Selection==4)


input = inputdlg({'Enter length of link 1 (m)','Enter length of link 2 (m)','Enter length of link 3 (m)','Enter length of link 4 (m)','Enter division ratio (r:1) of point on coupler'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});
l1=str2num(input{1});
l2=str2num(input{2});
l3=str2num(input{3});
l4=str2num(input{4});
r=str2num(input{5});


maxl=max([l1,l2,l3,l4]);
maxl=maxl*3;


limit1=0;
limit2=2*pi;

if(min([l1,l2,l3,l4])==l1&&l2+l3+l4-max([l1,l2,l3,l4])>=max([l1,l2,l3,l4])+l1)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DOUBLECRANK%%%%%

for(t2=0:0.01*pi:4*pi)

    %%POSITION ANALYSIS%%
[x1,y1]=circcirc(l2*cos(t2),l2*sin(t2),l3,l1,0,l4);
   if(abs(x1(1,1))<0.001)
        x1(1,1)=0;
    end
    if(abs(x1(1,2))<0.001)
        x1(1,2)=0;
    end
    if(abs(y1(1,1))<0.001)
        y1(1,1)=0;
    end
    if(abs(y1(1,2))<0.001)
        y1(1,2)=0;
    end

 
  
 if((x1(1,1))~=0&&(y1(1,1))~=0)%%configuration
     x=x1(1,1);
     y=y1(1,1);
 else
    x=x1(1,2);
   y=y1(1,2);
 end
 

 t3=atan((y-l2*sin(t2))/(x-l2*cos(t2)));
  t4=atan((y)/(x-l1)); 
  %%POSITION ANALYSIS END%%
 
  %%%%%%%%%%%%%%%%
  if(0)  
  %%VELCOITY ANALYSIS%%
  va=l2*w2;
  va_theta=t2+pi/2;
  
  vab_slope=tan(t3+pi/2);
  vb_slope=tan(t4+pi/2);
  
  A = [(1) (-vab_slope);
      (1) (-vb_slope)];
  B=[(l2*sin(t2)-vab_slope*l2*cos(t2));
      (0-vb_slope*l1)];
  soln=A\B;

  if(l2==l4)
     soln(1,1)=va*sin(va_theta);
     soln(2,1)= va*cos(va_theta);
  end;
  
  %%VELOCITY ANALYSIS END%%
   %%Acceleration Analysis%%
  aat=-1*l2*a2;
  aan=-1*l2*w2*w2;
  aat_theta=va_theta;
  aan_theta=t2;  
  aa_net_x = aat*cos(aat_theta)+aan*cos(aan_theta);
   aa_net_y = aat*sin(aat_theta)+aan*sin(aan_theta);
   

  aban= -((va*cos(va_theta)-soln(1,1))^2 + (va*sin(va_theta)-soln(2,1))^2)/l3;  
  acn= (soln(1,1)^2 + soln(2,1)^2) / l4;
  aban_theta=t3;
  acn_theta=t4;
  
  %%acceleration Analysis done%%
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%PLOT%%%%%
subplot(1,2,1),
hold off;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([0,l2*cos(t2)],[0,l2*sin(t2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);

hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([l2*cos(t2),x],[l2*sin(t2),y],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);
hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([l1,x],[0,y],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);
pt=['(' num2str(x) ', ' num2str(y) ')' ];
text(x+0.1,y+0.1,pt,'fontsize',10);
hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([0,l1],[0,0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);


subplot(1,2,2),
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
hold on;
plot((r*x+l2*cos(t2))/(r+1),(r*y+l2*sin(t2))/(r+1),'-rs','LineWidth',1,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',1);
axis manual;
axis([-maxl maxl -maxl maxl]);

pause(0.001);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
%%Velocity PLOT%%
subplot(2,2,2),
hold off;
axis manual;
axis([-1 1 -1 1]);
quiver(0,0,va*cos(va_theta),va*sin(va_theta));
hold on;
axis manual;
axis([-1 1 -1 1]);
quiver(0,0,soln(2,1),soln(1,1));
axis manual;
axis([-1 1 -1 1]);
hold on;
quiver(va*cos(va_theta),va*sin(va_theta),soln(2,1)-va*cos(va_theta),soln(1,1)-va*sin(va_theta));
axis manual;
axis([-1 1 -1 1]);
hold on;
%%Acceleration Plot%%
 
subplot(2,2,3),
hold off;
axis manual;
axis([-3 3 -3 3]);
quiver(0,0,aa_net_x,aa_net_y);
hold on;
axis manual;
axis([-3 3 -3 3]);
%quiver(aa_net_x,aa_net_y,aban*cos(aban_theta),aban*sin(aban_theta));
hold on;
axis manual;
axis([-3 3 -3 3]);
%quiver(0,0,acn*cos(acn_theta),acn*sin(acn_theta));
hold on;
pause(0.1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
display('dbl');


else
    
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CrankRocker%%%%%%%%%%
 
t3p=NaN;
xp=NaN;
t2=0;

for(t2=0.01:0.01*pi:2*pi-0.01)

[x1,y1]=circcirc(l2*cos(t2),l2*sin(t2),l3,l1,0,l4);


 if(x1(1,2)==0||y1(1,2)==0)
     x=x1(1,1);
     y=y1(1,1);
 else
     x=x1(1,2);
     y=y1(1,2);
 end
   if(isnan(x1(1,1))==1&&isnan(x1(1,2))==1&&isnan(y1(1,1))==1&&isnan(y1(1,2))==1)
       break;
   end
end
display(t2);

t=t2;

for(t2=t-0.01*pi:-0.01*pi:-t+0.01*pi)

[x1,y1]=circcirc(l2*cos(t2),l2*sin(t2),l3,l1,0,l4);

 if(isnan(x1(1,1))==0)
     x=x1(1,1);
     y=y1(1,1);
     
 else
    x=x1(1,2);
    y=y1(1,2);
end

  
  t3=atan((y-l2*sin(t2))/(x-l2*cos(t2)));
  t4=atan((y)/(x-l1)); 
  %%POSITION ANALYSIS END%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
 
%%%%%PLOT%%%%%
%%%%SIMULATION%%%

subplot(1,2,1),
hold off;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([0,l2*cos(t2)],[0,l2*sin(t2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);
hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([l2*cos(t2),x],[l2*sin(t2),y],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);
hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([l1,x],[0,y],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);

hold on;
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
plot([0,l1],[0,0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3);

subplot(1,2,2),
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);
hold on;
plot((r*x+l2*cos(t2))/(r+1),(r*y+l2*sin(t2))/(r+1),'-rs','LineWidth',1,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',1);
axis manual;
axis equal;
axis([-maxl maxl -maxl maxl]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause(0.01);
%%%END SIMULATION%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Selection==2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

input = inputdlg({'Enter length of link 1 (m)','Enter length of link 2 (m)','Enter length of link 3 (m)','Enter length of link 4 (m)','Enter the input crank angle(degrees)'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});
a= str2num(input{2});
if isempty(a)
    a = 10;
end


b=str2num(input{3});
if isempty(b)
    b = 10;
end

c=str2num(input{4});
if isempty(c)
    c = 10;
end

d=str2num(input{1});
if isempty(d)
    d = 10;
end

t2=str2num(input{5});
if isempty(t2)
    t2 = 10;
end
t2=t2*2*pi/360;
input = inputdlg({'Enter magnitude of force on Link4 (m)','Enter distance of the force from ground,along link4 (m)','Enter inclination of force wrt positive x in degrees'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});


Force=str2num(input{1});
if isempty(Force)
    Force = 1;
end
distance=str2num(input{2});
if isempty(distance)
    distance = 1;
end
T=str2num(input{3});
if isempty(T)
    T = 1;
end
T=T*2*pi/360;

if(0)
a=0.40;
b=1;
c=0.75;
d=0.5;
t2=2.0943951;
Force=80;
T=2.61799388;
distance=0.35;
end

if(d>(a+b+c))
    msgbox('Incompatible link lengths!')    
else

    
    %Position Analysis
    z2 =a*exp(1i*t2);
    k1=d/a;
    k2=d/c;
    k3=(a^2+c^2+d^2-b^2)/(2*a*c);
    k4=d/b;
    k5=((c^2-d^2-a^2-b^2)/(2*a*b));
    A=cos(t2)-k1-(k2*cos(t2))+k3;
    B=(-2*sin(t2));
    C=k1-((k2+1)*cos(t2))+k3;
    D=cos(t2)-k1+(k4*cos(t2))+k5;
    E=(-2*sin(t2));
    F=k1+((k4-1)*cos(t2))+k5;
    if((B^2-(4*A*C))<0)
       
    else
        if((a==b&&b==c&&c==d)||(a==c&&b==d))
            t4=t2;
        else
    t4=2*atan((-B-sqrt(B^2-(4*A*C)))/(2*A));%Computing angle theta4
        end
    t3=2*atan((-E-sqrt(E^2-(4*D*F)))/(2*D));%Computing angle theta3
    
    z4=c*exp(1i*t4);
  

     %%%%%%%%Plotting Code%%%%%%%%%
     
    %Plotting Position Vector Configuration Diagram
    plot([0 real(z2)],[0 imag(z2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
    hold on;
    axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    plot([real(z2) (d+real(z4))],[imag(z2) imag(z4)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold on;
    plot([(d+real(z4)) d],[imag(z4) 0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    quiver(a+distance*cos(t4),distance*sin(t4),0.3*cos(T),0.3*sin(T));
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold off;
    title('Simulation') ;
    text(0,-0.1,'O2','fontsize',10);%naming joints
    text(real(z2),imag(z2),'A','fontsize',10);
    text(d+real(z4),imag(z4),'B','fontsize',10);
    text(d,-0.1,'O4','fontsize',10);
    
    Fx=(-1*(Force*cos(T)*distance*sin(t4)-Force*sin(T)*distance*cos(t4)))/(cos(t3)*c*sin(t4)-sin(t3)*c*cos(t4));
    Torque = Fx*cos(t3) * a * sin(t2)-Fx*sin(t3) * a * cos(t2);
    msgbox([num2str(Torque),'N/m']);
     
    end
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Selection==5)
    
input = inputdlg({'Enter length of link 1 (m)','Enter length of link 2 (m)','Enter length of link 3 (m)','Enter length of link 4 (m)','Enter the input crank angle in degrees.'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});
a= str2num(input{2});
if isempty(a)
    a = 10;
end


b=str2num(input{3});
if isempty(b)
    b = 10;
end

c=str2num(input{4});
if isempty(c)
    c = 10;
     
end

d=str2num(input{1});
if isempty(d)
    d = 10;
end

t2=str2num(input{5});
if isempty(t2)
    t2 = 10;
end

t2=t2*2*pi/360;

input = inputdlg({'Enter magnitude of Torque on Link1 in Nm','Enter distance of the output force from ground,along link4 (m)'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1'});


Torque=str2num(input{1});
if isempty(Torque)
    Torque = 1;
end
distance=str2num(input{2});
if isempty(distance)
    distance = 1;
end



if(0)
a=0.40;
b=1;
c=0.75;
d=0.5;
t2=2.0943951;
Torque=18.9;
T=2.61799388;
distance=0.35;
end

if(d>(a+b+c))
    msgbox('Incompatible link lengths!')    
else

    
    %Position Analysis
    z2 =a*exp(1i*t2);
    k1=d/a;
    k2=d/c;
    k3=(a^2+c^2+d^2-b^2)/(2*a*c);
    k4=d/b;
    k5=((c^2-d^2-a^2-b^2)/(2*a*b));
    A=cos(t2)-k1-(k2*cos(t2))+k3;
    B=(-2*sin(t2));
    C=k1-((k2+1)*cos(t2))+k3;
    D=cos(t2)-k1+(k4*cos(t2))+k5;
    E=(-2*sin(t2));
    F=k1+((k4-1)*cos(t2))+k5;
    if((B^2-(4*A*C))<0)
       
    else
        if((a==b&&b==c&&c==d)||(a==c&&b==d))
            t4=t2;
        else
    t4=2*atan((-B-sqrt(B^2-(4*A*C)))/(2*A));%Computing angle theta4
        end
    t3=2*atan((-E-sqrt(E^2-(4*D*F)))/(2*D));%Computing angle theta3
    
    z4=c*exp(1i*t4);
  T= 2.61799388;

     %%%%%%%%Plotting Code%%%%%%%%%
     
    %Plotting Position Vector Configuration Diagram
    plot([0 real(z2)],[0 imag(z2)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
    hold on;
    axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    plot([real(z2) (d+real(z4))],[imag(z2) imag(z4)],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold on;
    plot([(d+real(z4)) d],[imag(z4) 0],'-rs','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',3)
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    quiver(a+distance*cos(t4),distance*sin(t4),0.3*cos(T),0.3*sin(T));
      axis([-a-1 d+c+1 -(c+a) c+a])  ;
    axis manual;
    hold off;
    title('Simulation') ;
    text(0,-0.1,'O2','fontsize',10);%naming joints
    text(real(z2),imag(z2),'A','fontsize',10);
    text(d+real(z4),imag(z4),'B','fontsize',10);
    text(d,-0.1,'O4','fontsize',10);
    Fx=Torque/(cos(t3) * a * sin(t2)-sin(t3) * a * cos(t2));
   
   Force=(-1* Fx*(cos(t3)*c*sin(t4)-sin(t3)*c*cos(t4)))/(cos(T)*distance*sin(t4)-sin(T)*distance*cos(t4));
    
    msgbox(['Output Force: ',num2str(Force),'N']);
     
    end
end
end

if(Selection==3)


input = inputdlg({'enter length of ground link','enter lenght of input link','enter length of coupler link','enter length of output link','enter input angle','enter the mass of input link','enter the mass of coupler link','enter the mass of output link','enter the mass moment of inertia of link 2','enter the mass moment of inertia of link 3','enter the mass moment of inertia of link 4','enter angular acceleration of input link','enter angular velocity'},'Four bar mechanism analysis. M2-Group B',1,{'1','1','1','1','1','1','1','1','1','1','1','1','1'});

l1=str2num(input{1});
l2=str2num(input{2});
l3=str2num(input{3});
l4=str2num(input{4});
t12p=str2num(input{5});
m2=str2num(input{6});
m3=str2num(input{7});
m4=str2num(input{8});
I2=str2num(input{9});
I3=str2num(input{10});
I4=str2num(input{11});
alph12=str2num(input{12});
w12=str2num(input{13});
    
if(0)
l1=1;
 l2=.8;
 l3=.8;
 l4=.44;
 t12p=37;
 w12=30;
 alph12=3;
 m2=10;
 m4=4;
 m3=8;
 I2=.013;
 I3=.0135;
 I4=.026;
end

t12=t12p*pi/180;
r=sqrt((l2*cos(t12)-l1)^2+(l2*sin(t12))^2);
t_lamda=acos((l4^2+r^2-l3^2)/(2*l4*r));
t_mu=acos((l4^2+l3^2-r^2)/(2*l4*l3));
t_phi=acos((l2*cos(t12)-l1)/r);
t14=t_phi-t_lamda;
t13=t14-t_mu;
%velocity analysis
w13=w12*(l2*sin(t12-t14))/(l3*sin(t14-t13));
w14=w12*(l2*sin(t12-t13))/(l4*sin(t14-t13));
%acceleration analysis
alph13=(l2*w12^2*cos(t12-t14)+l2*alph12*sin(t12-t14)+l3*w13^2*cos(t14-t13)-l4*w14^2)/(l3*sin(t14-t13));
alph14=(l2*w12^2*cos(t12-t13)+l2*alph12*sin(t12-t14)-l4*w14^2*cos(t14-t13)+l3*w13^2)/(l4*sin(t14-t13));
aG2x=-(l2/2)*w12^2*cos(t12)-(l2/2)*alph12*sin(t12);
aG2y=-(l2/2)*w12^2*sin(t12)+(l2/2)*alph12*cos(t12);
aG3x=-l2*w12^2*cos(t12)-l2*alph12*sin(t12)-(l3/2)*w13^2*cos(t13)-(l3/2)*alph13*sin(t13);
aG3y=-l2*w12^2*sin(t12)+l2*alph12*cos(t12)-(l3/2)*w13^2*sin(t13)+(l3/2)*alph13*cos(t13);
aG4x=-(l4/2)*w14^2*cos(t14)-(l4/2)*alph14*sin(t14);
aG4y=-(l4/2)*w14^2*sin(t14)+(l4/2)*alph14*cos(t14);
r2x=l2*cos(t12);
r2y=l2*sin(t13);
r3x=l3*cos(t13);
r3y=l3*sin(t13);
r4x=l4*cos(t14);
r4y=l4*sin(t14);
A=[-1 0 1 0 0 0 0 0 0;0 -1 0 1 0 0 0 0 0;0 0 -r2y r2x 0 0 0 0 1;0 0 -1 0 1 0 0 0 0;0 0 0 -1 0 1 0 0 0;0 0 0 0 -r3y r3x 0 0 0;0 0 0 0 -1 0 1 0 0;0 0 0 0 0 -1 0 1 0;0 0 0 0 0 0 -r4y r4x 0];
B=[m2*aG2x;m2*aG2y;I2*alph12;m3*aG3x;m3*aG3y;I3*alph13;m4*aG4x;m4*aG4y;I4*alph14];
X=A\B;
F21x=X(1,1);
F21y=X(2,1);
F32x=X(3,1);
F32y=X(4,1);
F43x=X(5,1);
F43y=X(6,1);
F14x=X(7,1);
F14y=X(8,1);
Tau=X(9,1);
msgbox({'Force Values are in N' 'Torque is in Nm' 'F21x=' num2str(F21x)  'F21y=' num2str(F21y) 'F32x=' num2str(F32x) 'F32y=' num2str(F32y) 'F43x=' num2str(F43x) 'F43y=' num2str(F43y) 'F14x=' num2str(F14x) 'F14y=' num2str(F14y) 'Ts=' num2str(Tau)});
axis([-50 50 -50 50]);
title('inertia forces');
axis([-50 50 -50 50]);
x1=0;
y1=0;
x2=l2*cos(t12);
y2=l2*cos(t12);
x4=l1;
y4=0;
x3=l1+l4*cos(t14);
y3=l4*sin(t14); 
       
      subplot(1,2,1);
      title('Forces due to inertia');
      plot([x2 x3],[y2 y3],'r','linewidth',1);
       hold on;
        plot([x3 x4],[y3 y4],'g','linewidth',1);
        hold on;
        plot([x1 x4],[y1 y4],'y','linewidth',1);
       hold on;    
        plot([x1 x2],[y1 y2],'b','linewidth',1);
        hold on;
        Fi2x=-m2*aG2x;
        Fi2y=-m2*aG2y;
        Fi3x=-m3*aG3x;
        Fi3y=-m3*aG3y;
        Fi4x=-m4*aG4x;
        Fi4y=-m4*aG4y;
       x12=(x1+x2)/2;y12=(y1+y2)/2;
       x23=(x2+x3)/2;y23=(y2+y3)/2;
       x34=(x3+x4)/2;y34=(y3+y4)/2;
       quiver(x12,y12,Fi2x,Fi2y,.0001);
       hold on;
       quiver(x23,y23,Fi3x,Fi3y,.0001);
       hold on;
       quiver(x34,y34,Fi4x,Fi4y,.0001);      
       
       subplot(1,2,2);title('Force Vectors at Joints');
       plot([x2 x3],[y2 y3],'r','linewidth',1);
       hold on;
        plot([x3 x4],[y3 y4],'g','linewidth',1);
        hold on;
        plot([x1 x4],[y1 y4],'y','linewidth',1);
       hold on;    
        plot([x1 x2],[y1 y2],'b','linewidth',1);
        hold on;
        quiver(x1,y1,F14x,F14y,.0001);
        hold on;
        quiver(x1,y1,-F14x,-F14y,.0001);
         hold on;
         quiver(x2,y2,F32x,F32y,.0001);
         hold on;
         quiver(x2,y2,-F32x,-F32y,.0001);
          hold on;
          quiver(x3,y3,F43x,F43y,.0001);
          hold on;
          quiver(x3,y3,-F43x,-F43y,.0001);
          hold on;
          quiver(x4,y4,F14x,F14y,.0001);
          hold on;
          quiver(x4,y4,-F14x,-F14y,.0001);         
          
          
end

 display(Selection);   
end



