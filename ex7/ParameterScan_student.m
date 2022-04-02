clc
close all
clear

% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chemin d'acces au code compile
repertoire = './'; % './' on Linux, '' on Windows
executable = 'Exercice7_2022_student.exe'; % Nom de l'executable
mydirectory= pwd;

Nx=101;
Ny=3;
paramstr ='Nx'; % Nom du parametre a scanner
param    = Nx;  % Valeurs du parametre a scanner
nsimul   = numel(param);
input    = 'configuration.in.example';
request  = 'single_scan';

%% Variant to scan Nx and Ny together:

% Nx = (10:10:50)+1;
% Ny = Nx;
% paramstr ='N'; % Nom du parametre a scanner
% param=[Nx;Ny]; %Nx and Ny the grid points arrays to scan
% nsimul   = size(param,2);
% input    = 'configuration.in.example';
% request  = 'double_scan';

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(nsimul,1);

for ii = 1:nsimul
    
    delete_sim([mydirectory,'/sim']); %Delete files in /sim directory
    if strcmp(request,'double_scan')
        % Variant to scan Nx and Ny together:
        eval(sprintf('!%s%s %s %s=%.15g %s=%.15g' , repertoire, executable, input, [paramstr,'x'], param(1,ii), [paramstr,'y'], param(2,ii)));
        ntime = size(dir('sim'),1)-2;
        output{ii} = generateHfield(ntime,Nx(ii),Ny(ii));
    else
        eval(sprintf('!%s%s %s %s=%.15g', repertoire, executable, input, paramstr, param(ii)));
        ntime = size(dir('sim'),1)-2;
        output{ii} = generateHfield(ntime,Nx,Ny);
    end
    disp('Done.')
    
end

%% Figures %%
%%%%%%%%%%%%%

for ii=1:nsimul
    figure(ii)
    time = output{ii}.time;
    X    = output{ii}.X;
    Y    = output{ii}.Y;
    H    = output{ii}.H;
    make_xt_plot(time,X,Y,H)
end






function delete_sim(sim)
% Specify the folder where the files live.
myFolder = sim;
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.out'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now deleting %s\n', fullFileName);
    delete(fullFileName);
end
end


function out = generateHfield(ntime,Nx,Ny)
H=zeros(Ny,Nx,ntime);
time=zeros(ntime,1);
for ii=1:ntime
    fichier    = ['sim/output.',num2str(ii),'.out'];
    data_str   = importdata(fichier,' ',1);
    time(ii)   = str2double(data_str.textdata{1});
    data       = data_str.data;
    H(:,:,ii)  = reshape(data(:,3),Ny,Nx);
end
%Gathering data in a struct
out.time = time;
out.X = data(1:Ny:Nx*Ny,1);
out.Y = data(1:Ny,2);
out.H = H;

end

function make_xt_plot(time,X,Y,H)
%find y/2 to plot t vs x at y=y/2
Ny = numel(Y);
index_y=round(Ny/2);
Hcut=squeeze(H(index_y,:,:))';
contourf(X,time,Hcut);
xlabel('x [m]')
ylabel('t [s]')
title(['H(x,y=',num2str(Y(index_y)),') [m]'])
colorbar
end