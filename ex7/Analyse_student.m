clc
close all
clear

% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=100;
Ny=3;
ntime = size(dir('sim'),1)-2;

animate = 1;% animate=1 to plot a sequence of snapshots of the solution
H=zeros(Ny,Nx,ntime);
time=zeros(ntime,1);

 for ii=1:ntime
    fichier    = ['sim/output.',num2str(ii),'.out'];
    data_str   = importdata(fichier,' ',1);
    time(ii)   = str2double(data_str.textdata{1});
    data       = data_str.data;
    H(:,:,ii)  = reshape(data(:,3),Ny,Nx);

 end

X         = data(1:Ny:Nx*Ny,1);
Y         = data(1:Ny,2);

%% Figures %%
%%%%%%%%%%%%%

if animate
    figure
    for ii=1:ntime
        contourf(X,Y,H(:,:,ii),15,'LineStyle','None')
        xlabel('x [m]')
        ylabel('y [m]')
        title('H(x,y) [m]')
        colorbar
        axis equal
        axis([min(X) max(X) min(Y) max(Y)])
        pause(.01)
    end
end
%%
figure

index_y=round(Ny/2);
Hcut=squeeze(H(index_y,:,:))';
contourf(X,time,Hcut);
xlabel('x [m]')
ylabel('t [s]')
title(['H(x,y=',num2str(Y(index_y)),') [m]'])
colorbar


