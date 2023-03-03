%%Parameter Definition
% lw = 0.055; %Whisker Length
lb = 2.5/1000; %Bottom plate edge length/2

%Known Paremeters
% HeightofPusher=30;
HeightOfExtrusion=5;
MaxPixelsXhalf=1280/2;
MaxPixelsYhalf=720/2;
AOVyhalf=pi/180*16.26;
AOVxhalf=pi/180*26.95;
% Corrector=MaxPixelsXhalf/MaxPixelsYhalf*tan(AOVyhalf)/tan(AOVxhalf);
DistAway=26;
% numSprings=4;
% ZeroRowOffset=0;
PixelDensity=1/2*((DistAway-HeightOfExtrusion)/MaxPixelsXhalf*tan(AOVxhalf)+(DistAway-HeightOfExtrusion)/MaxPixelsYhalf*tan(AOVyhalf));
numSprings=4;

Theta=0*3.14159;
Dz=0/1000;
figure()
% p = 0.004;
% l = 0.003;
l = 0.001;
h = .0008*10^3;
n=10;
w=0.00080;
widthLow=.25
widthHigh=2
thicknessLow=.25
thicknessHigh=2
l = (0.00375-w)/2;
p = 0.00141+w;
newpoints = 500;
% insert equation to solve for min Dz detectable
minDz=.06*10^-3

[w,t] = meshgrid(linspace(widthLow,widthHigh,newpoints),linspace(thicknessLow,thicknessHigh,newpoints));
i=1;
for wi=widthLow:(widthHigh-widthLow)/(newpoints-1):widthHigh
    j=0;
    for hi=thicknessLow:(thicknessHigh-thicknessLow)/(newpoints-1):thicknessHigh
        w1=wi/1000;
        h1=hi/1000;
        j=j+1;
        p = 0.00141+w1;
        Phi=0;
        [dbx,dby]=get_plate_disp2(0,0, minDz);
        [Mx, My, Fz] = get_momentPameter(0, 0, dbx, dby, minDz, lb,numSprings,h1,w1,n,l,p);
        BDmatrix(j,i)=Fz*1000;
        Fz;
    end
i=i+1;
end
%BDmatrixq = interp2('XData',linspace(widthLow,widthHigh,newpoints ),'YData',linspace(thicknessLow,thicknessHigh,newpoints) ,BDmatrix,w,t,'cubic');

%[c,h]=image('XData',linspace(widthLow,widthHigh,newpoints ),'YData',linspace(thicknessLow,thicknessHigh,newpoints),'CData',-BDmatrix);
%[C,h]=contourf(w,t, log(-BDmatrix),500)
x=linspace(widthLow,widthHigh,newpoints);
y=linspace(thicknessLow,thicknessHigh,newpoints);

imagesc(x,y,-transpose(BDmatrix))
hold on
lvls=[.001,.01]*1000;
[YY,h2]=contour(w,t, -transpose(BDmatrix),[4.0, 100],'LineColor',[.9 .9 .9],'LineStyle',['-'])
%clabel(C,h,'manual','FontSize',30,'Color','k')
hold on
[XX,h1]=contour(w,t, -transpose(BDmatrix),[1.05, 100],'LineColor','k','LineStyle',['--'])
%clabel(C,h,'manual','FontSize',30,'Color','k')
hold on
h1.LineWidth = 4;
h2.LineWidth = 4;
plot(.75, .75, 'p', 'LineWidth', 2, 'MarkerSize', 30,'MarkerFaceColor','r','MarkerEdgeColor','r')
%set(h,'LineColor','none')
%cbar = colorbar('Yscale', 'log');
colormap('parula')
a=colorbar
%a.Label.String = 'Minimum Detectable Force (N)'
%a.Title.String='Minimum Detectable Force (N)';
% ylabel(a,{'Minimum Detectable' ; 'Force (mN)'},'FontSize',30,'Rotation',270);
% pos = get(a,'Position');
% a.Label.Position=[8 2]
set(gca,'ColorScale','log')
hold on
% xlabel('Thickness (t) [mm]','FontSize',30)
% ylabel('Width (w) [mm]','FontSize',30)
ax = gca;
ax.FontSize = 36;
set(gca,'YDir','normal')

axis square

i=1;
figure()
for si=widthLow:(widthHigh-widthLow)/(newpoints-1):widthHigh
    j=0;
    for ri=.25:(20-.25)/(newpoints-1):20
        j=j+1;
        p = 0.00141+si/1000;
        Phi=0;
        %PixelDensity=1/2*((DistAway-ri)/MaxPixelsXhalf*tan(AOVxhalf)+(DistAway-ri)/MaxPixelsYhalf*tan(AOVyhalf));
        pi_w=atan(PixelDensity*1/ri);
        [theta_x, theta_y] = get_cart_angle(pi_w, 0);
        [dbx,dby]=get_plate_disp2(theta_x,theta_y, lb);
        [Mx, My, Fz] = get_momentPameter(theta_x, theta_y, dbx, dby, 0, lb,numSprings,si/1000,.75*10^-3,n,l,p);
        BDmatrix2(i,j)=sqrt(Mx^2+My^2);
        
        
    end
i=i+1;
end
x2=linspace(widthLow,widthHigh,newpoints);
y2=linspace(.25,20,newpoints);
[s,r] = meshgrid(x2,y2);
imagesc(x2,y2,transpose(BDmatrix2))
hold on
lvls=[.5];
[C,h]=contour(s,r, transpose(BDmatrix2),lvls,'LineColor',[.75 .75 .75])
clabel(C,h,'manual','FontSize',30,'Color','k')
hold on
h.LineWidth = 4;
colormap('parula')
a2=colorbar
% ylabel(a2,{'Minimum Detectable' ; 'Moment [N-mm]'},'FontSize',30,'Rotation',270);
% %a2.Label.String = 'Minimum Detectable' ;'Force (N)'
% pos = get(a2,'Position');
% a2.Label.Position=[8 .05]
set(gca,'ColorScale','log')
hold on
% xlabel('Thickness (t) [mm]','FontSize',30)
% ylabel('(r_h) [mm]','FontSize',30)
ax = gca;
ax.FontSize = 36;
plot(.75, 7.5, 'p', 'LineWidth', 2, 'MarkerSize', 30,'MarkerFaceColor','r','MarkerEdgeColor','r')
set(gca,'YDir','normal')
axis square

function [theta_x, theta_y] = get_cart_angle(pi_w, theta_w)
%     theta_y=asin((sin(pi_w)*sin(theta_w))/sqrt((sin(pi_w)*sin(theta_w))^2+cos(pi_w)^2))
%     theta_x=asin((sin(pi_w)*cos(theta_w))/sqrt((sin(pi_w)*cos(theta_w))^2+cos(pi_w)^2))
    theta_x = atan(sin(pi_w)*cos(theta_w)/cos(pi_w));
    theta_y = atan(sin(pi_w)*sin(theta_w)/cos(pi_w));
end
function [db1, db2] = get_plate_disp2(pi_w, theta_w, lb)
    W = [0 0 -cos(theta_w); 0 0 -sin(theta_w); cos(theta_w) sin(theta_w) 0];
    R = eye(3) + sin(pi_w)*W + 2*sin(pi_w/2)^2*W*W;
    x1 = [lb, 0, 0];
    y1 = [0, lb, 0];
    D1 = R*x1';
    D2 = R*y1';
    db2 = D1(3,1);
    db1 = D2(3,1);
end
function [db1, db2] = get_plate_disp(theta_x, theta_y, lb)
    db2=lb*tan(theta_x);
    db1=lb*tan(theta_y);
%     W = [0 0 -cos(theta_w); 0 0 -sin(theta_w); cos(theta_w) sin(theta_w) 0];
%     R = eye(3) + sin(pi_w)*W + 2*sin(pi_w/2)^2*W*W;
%     x1 = [lb, 0, 0];
%     y1 = [0, lb, 0];
%     D1 = R*x1';
%     D2 = R*y1';
%     db1 = D1(3,1);
%     db2 = D2(3,1);

end
function [Mx, My, Fz] = get_moment(theta_x, theta_y, dbx, dby, dTot, lb,numSprings)
%     %Define spring variables
%     p = 0.0003;
%     t = 0;
%     l = 0.0005;
%     w = 1.15*10^-4;
%     h = 0.0001;
%     n = 6;
%     E = 2*10^11; %% Young's modulus
%     v = 0.35; 
%     Gu = 77.2*10^9
%     %E/2/(1+v);
%     Ia= w*h*h*h/12;
%     Ib=Ia
% %     Ib = w*w*w*h/12;
%     Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
%     Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    %Define spring variables
    p = 0.0026;
    t = 0;
    l = 0.003;
    w = 0.00075;
    h = 0.00075;
    n = 10;
    E = 3.5*10^9; %% Young's modulus
    v = 0.37; 
    Gu = E/2/(1+v);
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    Ib=Ia;
    if w>h
        Jb=Ja;
    else
        Ja=Jb;
    end
    

    %%%Define Acquiring Values%%%
    Out=zeros(numSprings,3);

    %%%Equations for C Matrix%%%
    C = zeros(3);
    C(1,1) = 1/E/Ia*p*(2+n)+1/Gu/Jb*(n*(n-1)*t+2*n*l);
    C(1,2) = 0;
    C(2,1) = 0;
    C(1,3) = -1/E/Ia*p*p*(2+2*n+n*n/2)-1/Gu/Jb*p*(1/6*n*(4*n+7)*(n-1)*t+n*(n+2)*l);
    C(3,1) = C(1,3);
    C(2,2) = 1/Gu/Ja*p*(2+n)+1/E/Ib*((2*n+1)*l+n*(n-1)*t);
    C(2,3) = -1/Gu/Ja*p*t*n/2+1/E/Ib*(-n*(n-1)*t*t/2-n*l*t);
    C(3,2) = C(2,3);
    C(3,3) = 1/E/Ia*p*p*p*(1/3*(2+n)^3)+1/Gu/Ja*p*(1/6*n*(2*n-1)*(n-1)*t*t+n*(n-1)*l*t+n*l*l)+1/E/Ib*(1/6*n*n*(n-1)^2*t*t*t+2/3*n*l*l*l+n*(n-1)*l*l*t+(2/3*n*n*n-n*n+1/3*n)*l*t*t)+1/Gu/Jb*p*p*(1/2*n*(n-1)*(n*n+3*n+3)*t+(2/3*n*n*n+2*n*n+7/3*n)*l);
    
    %%Calculate Mx, My
    dSpring=[-dTot+dbx; -dTot+dby;-dTot-dbx;-dTot-dby];
    C1 = inv(C);
    for i=1:numSprings
        if rem(i,2)==0
            dis = [theta_x; theta_y; dSpring(i)];
        else
            dis = [theta_y; theta_x; dSpring(i)];
        end
        if i/2>1
            dis=dis.*[-1;-1;1];
        else
            dis;
        end
        Out(i,:)=C1*dis;
        %C2 = [C1(1,2) C1(1,1) C1(1,3) 0; C1(2,2) C1(2,1) C1(2,3) 0; C1(3,2) C1(3,1) C1(3,3) 0; C1(1,1) C1(1,2) 0 C1(1,3); C1(2,1) C1(2,2) 0 C1(2,3); C1(3,1) C1(3,2) 0 C1(3,3)];
    end
    Fz=sum(Out(:,3));
    Mx=Out(1,1)+Out(2,2)-Out(3,1)-Out(4,2)+lb*Out(1,3)-lb*Out(3,3);
    My=Out(2,1)+Out(1,2)-Out(4,1)-Out(3,2)+lb*Out(2,3)-lb*Out(4,3);
    Mx=Mx*1000;
    My=My*1000;
end

function [Mx, My, Fz] = get_momentPameter(theta_x, theta_y, dbx, dby, dTot, lb,numSprings,h,w,n,l,p)
%     %Define spring variables
%     p = 0.0003;
%     t = 0;
%     l = 0.0005;
%     w = 1.15*10^-4;
%     h = 0.0001;
%     n = 6;
%     E = 2*10^11; %% Young's modulus
%     v = 0.35; 
%     Gu = 77.2*10^9
%     %E/2/(1+v);
%     Ia= w*h*h*h/12;
%     Ib=Ia
% %     Ib = w*w*w*h/12;
%     Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
%     Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    %Define spring variables
    %p = 0.0015;
    t = 0;
    %l = 0.002;
    %w = 0.00075;
    %h = 0.0025;
    %n = 12;
    E = 2.1*10^9; %% Young's modulus
    v = 0.35; 
    Gu = E/2/(1+v);
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    Ib=Ia;
    if w>h
        Jb=Ja;
    else
        Ja=Jb;
    end
    
  

    %%%Define Acquiring Values%%%
    Out=zeros(numSprings,3);

    %%%Equations for C Matrix%%%
    C = zeros(3);
    C(1,1) = 1/E/Ia*p*(2+n)+1/Gu/Jb*(n*(n-1)*t+2*n*l);
    C(1,2) = 0;
    C(2,1) = 0;
    C(1,3) = -1/E/Ia*p*p*(2+2*n+n*n/2)-1/Gu/Jb*p*(1/6*n*(4*n+7)*(n-1)*t+n*(n+2)*l);
    C(3,1) = C(1,3);
    C(2,2) = 1/Gu/Ja*p*(2+n)+1/E/Ib*((2*n+1)*l+n*(n-1)*t);
    C(2,3) = -1/Gu/Ja*p*t*n/2+1/E/Ib*(-n*(n-1)*t*t/2-n*l*t);
    C(3,2) = C(2,3);
    C(3,3) = 1/E/Ia*p*p*p*(1/3*(2+n)^3)+1/Gu/Ja*p*(1/6*n*(2*n-1)*(n-1)*t*t+n*(n-1)*l*t+n*l*l)+1/E/Ib*(1/6*n*n*(n-1)^2*t*t*t+2/3*n*l*l*l+n*(n-1)*l*l*t+(2/3*n*n*n-n*n+1/3*n)*l*t*t)+1/Gu/Jb*p*p*(1/2*n*(n-1)*(n*n+3*n+3)*t+(2/3*n*n*n+2*n*n+7/3*n)*l);
    
    %%Calculate Mx, My
    dSpring=[-dTot+dbx; -dTot+dby;-dTot-dbx;-dTot-dby];
    C1 = inv(C);
    for i=1:numSprings
        if rem(i,2)==0
            dis = [theta_x; theta_y; dSpring(i)];
        else
            dis = [theta_y; theta_x; dSpring(i)];
        end
        if i/2>1
            dis=dis.*[-1;-1;1];
        else
            dis;
        end
        Out(i,:)=C1*dis;
        %C2 = [C1(1,2) C1(1,1) C1(1,3) 0; C1(2,2) C1(2,1) C1(2,3) 0; C1(3,2) C1(3,1) C1(3,3) 0; C1(1,1) C1(1,2) 0 C1(1,3); C1(2,1) C1(2,2) 0 C1(2,3); C1(3,1) C1(3,2) 0 C1(3,3)];
    end
    Fz=sum(Out(:,3));
    Mx=Out(1,1)+Out(2,2)-Out(3,1)-Out(4,2)+lb*Out(1,3)-lb*Out(3,3);
    My=Out(2,1)+Out(1,2)-Out(4,1)-Out(3,2)+lb*Out(2,3)-lb*Out(4,3);
    Mx=Mx*1000;
    My=My*1000;
end
