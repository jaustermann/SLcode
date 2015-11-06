% Code to solve the sea level equation. Follows class notes of Jerry
% Mitrovica (lecture 17 of sea level notes)
% Code solves the elastic case neglecting rotational effects, moving
% shorelines, viscoelastic behaviour, gravitational effects from ice sheets
% J. Austermann 2012

% Specify maximum degree to which spherical transformations should be done
maxdeg = 256;

% parameters
rho_ice = 920;
rho_water = 1000;

tic

% ice masks / contours
% save('disk3_icehist','ice_height','t');
load ice5g_icehist
% load disk3_icehist_10
timestep = t;
% time is from -21 to 0

% for k = 1:11
%     ice_height_k{k} = ice_height{11};
% end
% 
% ice_height = ice_height_k;

% time normalization because spoles of maxwell are outputted in 1/1000yrs
timestep = timestep/1000;

LatKonst = double(ice_lat' + 90);
LonKonst = double(ice_lon');
ice_height = ice_height_ice5g;
% for disk runs
% LatKonst = 180:-1:0;
% LonKonst = 0:1:360;

ice_lm = cell(length(timestep));
dice_lm = cell(length(timestep));

for i = 1:length(timestep)
    ice_lm{i} =  spa2sph(ice_height{i},maxdeg,LonKonst,LatKonst);
    
    if i == 1
        ice_diff = zeros(size(ice_height{1}));
    else
        ice_diff = ice_height{i} - ice_height{i-1};
    end
    dice_lm{i} = spa2sph(ice_diff,maxdeg,LonKonst,LatKonst);

end

% prepare love numbers in suitable format and calculate T_lm and E_lm 
load ReadSaveLN/LN_l90_VM2_Jerry
% load ReadSaveLN/LN_l90_VM2L3
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);

% E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);


% Calculate viscous component
V_lm_h = cell(length(timestep));
V_lm_k = cell(length(timestep));
R_lm = cell(length(timestep));
G_lm = cell(length(timestep));
R = cell(length(timestep));
G = cell(length(timestep));

V_lm_h{1} = zeros(size(T_lm));
V_lm_k{1} = zeros(size(T_lm));

for j = 2:length(timestep)
    
    visc_h = zeros(j-1,length(h_lm));
    visc_k = zeros(j-1,length(h_lm));
    
    for n = 1:j-1
        
        beta_h = zeros(maxdeg, 1);
        beta_k = zeros(maxdeg, 1);
        for lm = 1:maxdeg
            k = mode_found(lm);
            beta_h(lm) = sum(h_amp(lm,1:k)./spoles(lm,1:k).* ...
            (1 - exp(- spoles(lm,1:k) * (timestep(j) - timestep(n)))));
            beta_k(lm) = sum(k_amp(lm,1:k)./spoles(lm,1:k).* ...
            (1 - exp(- spoles(lm,1:k) * (timestep(j) - timestep(n)))));
        end
        beta_lm_h = love_lm(beta_h, maxdeg);
        visc_h(n, 1:length(beta_lm_h)) = beta_lm_h .* (rho_ice * dice_lm{n}); 
        
        beta_lm_k = love_lm(beta_k, maxdeg);
        visc_k(n, 1:length(beta_lm_k)) = beta_lm_k .* (rho_ice * dice_lm{n}); 
    end
    
    V_lm_h{j} = sum(visc_h,1);
    V_lm_k{j} = sum(visc_k,1);
    
    disp(['First loop done by ' num2str(j/length(timestep)*100) '%'])
end


for i = 1:length(timestep)
    R_lm{i} = T_lm.*(h_lm.*(rho_ice*(ice_lm{i}-ice_lm{1})) + V_lm_h{i});
    R_lm{i}(1) = 0; %nodcshift
    G_lm{i} = T_lm.*((1+k_lm).*(rho_ice*(ice_lm{i}-ice_lm{1})) + V_lm_k{i});
    G_lm{i}(1) = 0; %nodcshift
    % R_lm{i} = T_lm.*(E_lm.*(rho_ice*(ice_lm{i}-ice_lm{1})));
    % R{i} = sph2spa(R_lm{i},maxdeg,LonKonst,LatKonst');
    G{i} = sph2spa(G_lm{i},maxdeg,LonKonst,LatKonst');
    disp(['Second loop done by ' num2str(i/length(timestep)*100) '%'])
end

save('ice5g_res_step_LNJerry_G','R_lm','G_lm','R','G');


toc

%% plot result
% figure
% pcolor(LonKonst,LatKonst,R{2})
% shading flat
% colorbar

% figure
% subplot(2,1,1)
% title(['At time' num2str(t(1))])
% plot(180-LatKonst,R{1}(:,1))
% ylabel('Radial Displacement [m]')
% xlabel('Colatitude')
% subplot(2,1,2)
% plot(180-LatKonst,G{1}(:,1) * 9.81)
% ylabel('Grav. Potential [m^2/s^2]')
% xlabel('Colatitude')
% 
% figure
% subplot(2,1,1)
% title(['At time' num2str(t(11))])
% plot(180-LatKonst,R{11}(:,1))
% ylabel('Radial Displacement [m]')
% xlabel('Colatitude')
% subplot(2,1,2)
% plot(180-LatKonst,G{11}(:,1) * 9.81)
% ylabel('Grav. Potential [m^2/s^2]')
% xlabel('Colatitude')

% figure
% subplot(2,1,1)
% val = [];
% for i=1:length(timestep)
%     val(i) = R{i}(1,1);
% end
% plot(t,val)
% 
% subplot(2,1,2)
% val = [];
% for i=1:length(timestep)
%     val(i) = G{i}(1,1)  * 9.81;
% end
% plot(t,val)


% %Compare Konstantin
% 
% ind = find(t == 5000);
% figure
% rate_R = (R{ind+1}(:,1) - R{ind}(:,1))/(t(ind+1) - t(ind));
% 
% plot(180-LatKonst(1:71), rate_R(1:71)*1000);

%% spatial mesh
% i = 70;
% [LON,LAT] = meshgrid(LonKonst,LatKonst-90);
% figure
% subplot(2,1,1)
% pcolor(LON,LAT,ice_height{i} - ice_height{i-1})
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(LON,LAT,R{i})
% shading flat
% colorbar

%% extract point in hudson bay

% lon_HB = 270;
% lat_HB = 90-30;
% 
% for i = 1:length(timestep)
%     point_HB(i) = interp2(LON_LN,LAT_LN,R{i},lon_HB,lat_HB);
% end
% 
% figure
% plot(timestep,point_HB)


