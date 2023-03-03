load('map_E3D_tk_r_v1.mat')
[Dir,Con,Mag]=size(MZ);
Thetas=atan2(y_pt,z_pt);
Dist=(y_pt.^2+z_pt.^2).^(1/2);
Phi=atan(Dist./x_pt)*180/3.14159;
figure()
marker=['o','d','<'];
MSUM=sqrt(MZ.^2+MY.^2);
MD=atan2(MY,MZ);
graphColor=hsv(8);
% for i=1:Dir
%     for j=1:Con
%         for k=1:Mag
%             scatter3(MSUM(i,j,k)*1000,FX(i,j,k), MD(i,j,k),'o','MarkerFaceColor',graphColor(j,:),'MarkerEdgeColor','k')
%             hold on
%         end
%     end
% end
% xlabel('Moment (N-mm)','FontSize',16)
% ylabel('Force (N)/Moment (N-mm)','FontSize',16)
% zlabel('Direction Moment (deg)','FontSize',16)
% %axis([0 .003 -.01 .01])
hold off
figure()
marker=['o','^','^','^','^','^']
for i=[1]
    for j=1:Con
        for k=1:Mag
            plot((MSUM(i,j,k)*1000),-FX(i,j,k)*1000./(MSUM(i,j,k)*1000),marker(i),'MarkerFaceColor',graphColor(j,:),'MarkerEdgeColor','k','MarkerSize',8)
            hold on
        end
    end
end
ylabel('Fz (N )','FontSize',24)
xlabel('MB (N-mm)','FontSize',24)
patch([0 35 35 0], [0 0 .05 .05], [0.8 0.8 0.8])
% ylim([0 180])
% xlim([0 35])
set(gca,'FontSize',24)
%axis([0 .003 -.01 .01])
hold off
figure()

for i=[1]
    for j=1:Con
        for k=1:Mag
            plot(k,k/(-FX(i,j,k)),marker(i),'MarkerFaceColor',graphColor(j,:),'MarkerEdgeColor','k','MarkerSize',8)
            hold on
        end
    end
end
xlabel('translation [mm]','FontSize',24)
ylabel('Force [mN]','FontSize',24)
patch([0 180 180 0], [0 0 3 3], [0.8 0.8 0.8])
hold on
patch([0 .05 .05 0], [0 0 6000 6000], [0.8 0.8 0.8])
hold on
set(gca,'FontSize',24)
% xlim([1 180])
% ylim([1 6000])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
yline(3)
hold off

figure()
for i=[1,3,5,7]
    for j=1:Con-1
        for k=1:Mag
            plot(Phi(i,j,k),MX(i,j,k),marker(i),'MarkerFaceColor',graphColor(i,:),'MarkerEdgeColor','k')
            %plot(Phi(i,j,k),MZ(i,j,k),marker(j),'MarkerFaceColor',graphColor(i,:),'MarkerEdgeColor','k')
            hold on
        end
    end
end
xlabel('Phi Deflection (rad)')
ylabel('Moment (N-m)')
hold off

figure()
for i=1
    for j=1:Con
        for k=1:Mag
            scatter3(MZ(i,j,k),MY(i,j,k),MZ(i,j,k), 15,graphColor(i,:),'filled',marker(j))
            %plot(Phi(i,j,k),MZ(i,j,k),marker(j),'MarkerFaceColor',graphColor(i,:),'MarkerEdgeColor','k')
            hold on
        end
    end
end
xlabel('Moment Y (N-m)')
ylabel('Moment X (N-m)')
zlabel('Z Force (N)')
hold off
