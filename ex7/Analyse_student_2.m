%clc
%close all
%clear

% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=100; % MUST match the one used in the C++ code
Ny=3;  % MUST match the one used in the C++ code

ntime_all = size(dir('sim'),1)-2; % all time steps simulated
ntime_start=1; % start index of time for reading output files
ntime_end=ntime_all;%ntime_all; % end index of time for reading output files
istride=1; % reads only every istride time steps (>1 useful for large and long simulations)
% NB: if you want to read just one time ("snapshot"), put
% ntime_start=ntime_end

index_y=round(Ny/2); % cut at y=L_y/2 (change this according to your needs)

% different types of plots
% 'xyt': 2D contours f(x,y,t_ii) at timesteps ii=nstime_start..ntime_end
% 'xt': 1D f(x,y=ycut,t_ii) and 2D contours f(x,y=ycut,t)

nsel_plottype='xyt' % 

animate = 1; %0: plot only the ntime_end ; 1: animate from ntime_start to ntime_endS
nsteps = ntime_end-ntime_start+1;
if animate
    istart=1;
else
    istart=nsteps;
end
switch nsel_plottype
    case 'xyt'
        H=zeros(Ny,Nx,nsteps);
    case 'xt'
        H=zeros(Nx,nsteps);
end
time=zeros(nsteps,1);

 for ii=ntime_start:istride:ntime_end
    it=ii-ntime_start+1; % index of arrays
    fichier    = ['sim/output.',num2str(ii),'.out'];
    data_str   = importdata(fichier,' ',1); % reads the file and puts it into a data structure
    time(it)   = str2double(data_str.textdata{1});
    data       = data_str.data;
    switch nsel_plottype
        case 'xyt'
            H(:,:,it) = reshape(data(:,3),Ny,Nx);
        case 'xt'
            Hxy = reshape(data(:,3),Ny,Nx);
            H(:,it)= Hxy(index_y,:);
    end
 end

X         = data(1:Ny:Nx*Ny,1);
Y         = data(1:Ny,2);
dx        = X(2)-X(1);
dy        = Y(2)-Y(1);

%% Figures %%
%%%%%%%%%%%%%
lw=1; fs=16;
switch nsel_plottype

    case 'xyt'
        
        figure
        for ii=istart:nsteps
            contourf(X,Y,H(:,:,ii),15,'LineStyle','None')
            set(gca,'fontsize',fs)
            xlabel('x [m]')
            ylabel('y [m]')
            title('f(x,y)')
            colorbar
            axis equal
            axis([min(X) max(X) min(Y) max(Y)])
            pause(.01)
        end
        if nsteps>1  
            figure
            Hcut=squeeze(H(index_y,:,:))';
            contourf(X,time,Hcut);
            xlabel('x [m]')
            ylabel('t [s]')
            title(['f(x,y=',num2str(Y(index_y)),') [m]'])
            colorbar
        end

    case 'xt'
     
        figure
        Hmin=min(min(H)); Hmax=max(max(H));
        for ii=istart:nsteps
            plot(X,H(:,ii),'b-','linewidth',lw)
            set(gca,'fontsize',fs)
            xlabel('x [m]')
            ylabel('f')
            title(['f(x,y=',num2str(Y(index_y),3),')'])
            axis([min(X) max(X) Hmin Hmax])
            pause(.05)
        end
        if nsteps>1  
            figure
            contourf(X,time,H');
            xlabel('x [m]')
            ylabel('t [s]')
            title(['f(x,y=',num2str(Y(index_y)),')'])
            colorbar
        end
end
        
        

