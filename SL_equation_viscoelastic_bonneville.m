% Code to solve the elastic sea level equation following 
% Kendall et al., 2005 and Austermann et al., 2015

% J. Austermann 2015

% add paths when run for the first time.
% addpath SLFunctions
% addpath SLFunctions/mtimesx
% addpath '/Users/jackyaustermann/Documents/MATLAB/m_map'

%% Parameters & Input 
% Specify maximum degree to which spherical transformations should be done
maxdeg = 512;

% Some options to choose from
include_lakes = 'n'; % choose between y (for yes) and n (for no)
include_rotation = 'n'; % choose between y (for yes) and n (for no)
include_ice_check = 'n'; % choose between y (for yes) and n (for no)
% include_100year_timestep = 'n'; % choose between y (for yes) and n (for no)

% parameters
rho_ice = 920;
rho_water = 1000;
rho_sed = 2300;
g = 9.80665;


% The following steps help speed up the calculations
% Set up Gauss Legendre grid onto which to interpolate all grids
N = maxdeg; 
[x,w] = GaussQuad(N);
x_GL = acos(x)*180/pi - 90;
lon_GL = linspace(0,360,2*N+1);
lon_GL = lon_GL(1:end-1);

colat = 90 - x_GL;
lon = lon_GL;

[lon_out,lat_out] = meshgrid(lon_GL,x_GL);

% Precompute legendre polynomials
[P_lm_spa2sph, P_lm_sph2spa] = get_Legendre(x_GL,maxdeg);
% 
% 
% --------------------------------
% ICE
% --------------------------------

% load WAIS 
% load ice5g_griddata
% 
% if include_100year_timestep == 'n'
%     ice_time = [ice_time(1); ice_time(3:end)];
%     ice5g_grid = [ice5g_grid(1,:,:); ice5g_grid(3:end,:,:)];
% end

ice = zeros(length(x_GL),length(lon_GL),length(ice_time_new));
% ice_time_new = zeros(size(ice_time));
% 
% ice_lat = [90; ice_lat; -90];
% 
% for i = 1:length(ice_time)
% 
%     % already on Gauss Legendre grid
%     ice_nointerp = squeeze(ice5g_grid(i,:,:));
%     
%     % add rows at top and bottom
%     ice_extended = [zeros(1,length(ice_long)); ice_nointerp; ...
%         ice_nointerp(1,end)*ones(1,length(ice_long))];
% 
%     % interpolate ice masks on common grid
%     ice_interp = interp2(ice_long,ice_lat,ice_extended,lon_out, lat_out);
% 
%     % patch in zeros
%     % ice_interp(isnan(ice_interp) == 1) = 0;
%     
%     ice(:,:,length(ice_time)-i+1) = ice_interp;
%     ice_time_new(i) = ice_time(length(ice_time)-i+1);
% end

%r = find(ice_time_new==18);
%ice_time_new = ice_time_new(r:end);

ice_time_new = [0 1];


lakeB = zeros(length(x_GL),length(lon_GL),length(ice_time_new));

r_lat = find((38 < x_GL) & (x_GL < 42));
r_lon = find((360-114 < lon_GL) & (lon_GL < 360-112)); %111.5 - 114.5 lon

avg_depth = 350; %m

for i = 1:length(r_lat)
    for j = 1:length(r_lon)
        lakeB(r_lat(i),r_lon(j),1) = avg_depth;
    end
end

lake_00 = spa2sph(squeeze(lakeB(:,:,1)),...
    0,lon,colat,P_lm_spa2sph,w); %m3

area_00 = spa2sph(sign_01(-lakeB(:,:,1)+0.1),...
    0,lon,colat,P_lm_spa2sph,w);

lake_vol_calc = lake_00 * 4*pi*(6371e3)^2 * 1e-9; %km

ice = lakeB;
rho_ice = rho_water;


%area_lake = area_00 * 4*pi*(6371e3)^2 *1e-6; %km

%Vol - 10300 km^3
% depth 130m
% length 400km
% width 200km



% load ice5g_gl

% ice = cell(size(times));
% 
% for i = 1:length(times)
%     % interpolate ice masks on common grid
%     ice{i} = squeeze(ice5g_gl(:,1:end-1,i));
% end

% ice = single(ice5g_gl(:,1:end-1,:));
% clear ice5g_gl
% 
% ice_time_new = times;
% 
% 
% N = 256;%maxdeg; 
% [x,w] = GaussQuad(N);
% % x_GL = acos(x)*180/pi - 90;
% % lon_GL = linspace(0,360,2*N+1);
% % lon_GL = lon_GL(1:end-1);
% 
% x_GL = latitude;
% lon_GL = longitude(1:end-1);%lon_GL;
% 
% colat = 90 - x_GL;
% lon = lon_GL;
% 
% [lon_out,lat_out] = meshgrid(lon_GL,x_GL);
% 
% % Precompute legendre polynomials
% [P_lm_spa2sph, P_lm_sph2spa] = get_Legendre(x_GL,maxdeg);
% 
% P_lm_spa2sph = single(P_lm_spa2sph);
% P_lm_sph2spa = single(P_lm_sph2spa);
    

%del_ice = ice_j - ice_0; 


% --------------------------------
% DYNAMIC TOPOGRAPHY
% --------------------------------

%del_DT = zeros(size(lon_out));
DT = zeros(length(x_GL),length(lon_GL),length(ice_time_new));

% --------------------------------
% SEDIMENT
% --------------------------------

% sediment history 
sed = zeros(length(x_GL),length(lon_GL),length(ice_time_new));


% --------------------------------
% TOPOGRAPHY
% --------------------------------

% load preloaded etopo (including ice) as topo_orig, lon_topo, lat_topo
load topo_SL

% interpolate topography grid onto Gauss Legendre Grid
topo0 = interp2(lon_topo,lat_topo,topo_bed,lon_out, lat_out);% + ice{end};

% 
% load toposam.mat
% topo0 = topobr + ice(:,:,end);


%% Set up love number input

% prepare love numbers in suitable format and calculate T_lm and E_lm 
% to calculate the fluid case, switch h_el to h_fl, k_el to k_fl and same
% for tidal love numbers
%load SavedLN/prem.l90C.umVM2.lmVM2.mat
%load SavedLN/prem.l48C.ump2.lm7.mat
load SavedLN/prem.l21.ump1.lmVM2_coll_512.mat
h_lm = love_lm(h_fl, maxdeg);
k_lm = love_lm(k_fl, maxdeg);
% h_lm = love_lm(h_el, maxdeg);
% k_lm = love_lm(k_el, maxdeg);
h_lm_tide = love_lm(h_el_tide,maxdeg);
k_lm_tide = love_lm(k_el_tide,maxdeg);

E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);

E_lm_T = 1 + k_lm_tide - h_lm_tide;


% calculate time betas
beta_lm = cell(length(ice_time_new)-1);
% beta_lm = zeros(length(ice_time_new)-1,length(ice_time_new)-1,length(h_lm));

for t_it = 2:length(ice_time_new)
    
    for n = 2:t_it-1
        
        beta_o = zeros(maxdeg, 1);
        for lm = 1:maxdeg
            num_mod = mode_found(lm);
            beta_o(lm) = sum((k_amp(lm,1:num_mod) - h_amp(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n)))));
        end
        
%         mode_sum = (k_amp - h_amp)./spoles .* (1 - exp(- spoles ...
%                 * (-ice_time_new(t_it) + ice_time_new(n))));
%         beta = sum(mode_sum,2);
         beta_lm{t_it-1}(n-1,:) = single(love_lm(beta_o, maxdeg));

    end
end

%%
% calculate time betas
% beta_T_lm = zeros(length(ice_time_new)-1,length(ice_time_new)-1,length(h_lm));
% for t_it = 2:length(ice_time_new)
%     for n = 2:t_it-1
%         beta_T = zeros(maxdeg, 1);
%         for lm = 1:maxdeg
%             num_mod = mode_found(lm);
%             beta_T(lm) = sum((k_amp_tide(lm,1:num_mod) - h_amp_tide(lm,1:num_mod)) ...
%                 ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
%                 * (-ice_time_new(t_it) + ice_time_new(n)))));
%         end
%         beta_T_lm(t_it-1,n-1,:) = love_lm(beta_T, maxdeg);
%     end
% end
% 
% visc_k_L = zeros(length(ice_time_new)-1,length(ice_time_new)-1);
% visc_k_T = zeros(length(ice_time_new)-1,length(ice_time_new)-1);
% lm = 2;
% 
% for t_it = 2:length(ice_time_new)
%     for n = 1:t_it-1
%         num_mod = mode_found(lm);
%         visc_k_L(t_it-1,n-1,:) = sum((k_amp(lm,1:num_mod)) ...
%             ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
%             * (ice_time_new(t_it) - ice_time_new(n)))));
%         visc_k_T(t_it-1,n-1,:)) = sum((k_amp_tide(lm,1:num_mod)) ...
%             ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
%             * (ice_time_new(t_it) - ice_time_new(n)))));
%     end
% end
% 
% k_L = k_el(lm+1)+sum(visc_k_L,1);
% k_T = k_el(lm+1)+sum(visc_k_T,1);

%% Solve sea level equation (after Kendall 2005, Dalca 2013 & Austermann et al. 2015)
tic
k_max = 10;   % maximum number of iterations
epsilon = 10^-4; % convergence criterion

topo_it_max = 1;   % maximum number of iterations
max_topo_diff = 0.1; % convergence criterion

% 0 = before
% j = after

% set up initial topography and ocean function
% topo_initial = repmat(topo0,[1,1,topo_it_max+1]); % already includes ice, DT and sediments

% initial topography guess: topography is the same as present at every
% point in time; topography is a 3D vector; access topography at time x
% like this topo(:,:,x) [or for plotting squeeze(topo(:,:,x))]
% topo = repmat(topo0,[1,1,length(ice_time_new)]);

% set up initial topography and ocean function
topo_initial = zeros(length(x_GL),length(lon_GL),topo_it_max+1);
topo_initial(:,:,1) = topo0 - ice(:,:,end) + ice(:,:,1); % already includes ice, DT and sediments

% initial topography guess: topography is the same as present at every
% point in time; topography is a 3D vector; access topography at time x
% like this topo(:,:,x) [or for plotting squeeze(topo(:,:,x))]
topo = zeros(length(x_GL),length(lon_GL),length(ice_time_new));
for i = 2:length(ice_time_new)
    topo(:,:,i) = topo0 - ice(:,:,end) + ice(:,:,i);
end
% topo = repmat(topo0,[1,1,length(ice_time_new)]);


% initialize Lakes vector
Lakes = zeros(size(topo));
ice_corrected = zeros(size(topo));
sdelS_lm = zeros(length(ice_time_new),length(h_lm));
sdelL_lm = zeros(length(ice_time_new)-1,length(h_lm));

delS_lm_prev = zeros(1,length(h_lm));
TO_lm_prev = zeros(1,length(h_lm));
delL_lm_prev = zeros(1,length(h_lm));
delLa_lm_prev = zeros(1,length(h_lm));
deli_00_prev = 0;



% initial values for convergence
conv_topo = 'not converged yet';

for topo_it = 1:topo_it_max;
    
    topo_0 = topo_initial(:,:,topo_it);
    oc_0 = sign_01(topo_0);
    topo(:,:,1) = topo_0;
    ice_0 = ice(:,:,1);

    % expand ocean function into spherical harmonics
    oc0_lm = spa2sph(oc_0,maxdeg,lon,colat,P_lm_spa2sph,w);
    
    if include_ice_check == 'y'
        % check ice model for floating ice
        check1 = sign_01(-topo_0 + ice_0);
        check2 = sign_01(+topo_0 - ice_0) .* ...
            (sign_01(-ice_0*rho_ice - (topo_0 - ice_0)*rho_water));

        ice_0_corr = check1.*ice_0 + check2.*ice_0;
    else
        ice_0_corr = ice_0;
    end
    
    ice_corrected(:,:,1) = ice_0_corr;
    
    switch conv_topo

        case 'converged!'

        case 'not converged yet'

        for t_it = 2:length(ice_time_new) % loop over time
            
            % take topogaphy from last topo iteration as starting
            % topography
            topo_j = topo(:,:,t_it);
            oc_j = sign_01(topo_j);
            ocj_lm = spa2sph(oc_j,maxdeg,lon,colat,P_lm_spa2sph,w);
            ice_j = ice(:,:,t_it);

            % calculate change in sediments and decompose into spherical harmonics
            % set to zero as is but can be set as the change between the sed load. 
            del_sed = sed(:,:,t_it) - sed(:,:,1);
            delSed_lm = zeros(size(oc0_lm));% SWITCH IF USE SED spa2sph(del_sed,maxdeg,lon,colat,P_lm_spa2sph,w);
            
            del_DT = DT(:,:,t_it) - DT(:,:,1);
            delDT_lm = zeros(size(oc0_lm));% SWITCH IF USE SEDspa2sph(del_sed,maxdeg,lon,colat,P_lm_spa2sph,w);
            
            if include_ice_check == 'y'
                % check ice model for floating ice
                check1 = sign_01(-topo_j + ice_j);
                check2 = sign_01(+topo_j - ice_j) .* ...
                    (sign_01(-ice_j*rho_ice - (topo_j - ice_j)*rho_water));

                ice_j_corr = check1.*ice_j + check2.*ice_j;
            else
                ice_j_corr = ice_j;
            end
            
            ice_corrected(:,:,t_it) = ice_j_corr;
            
            del_ice_corrected = ice_corrected(:,:,t_it) - ice_0;
            deli_lm = spa2sph(del_ice_corrected,maxdeg,lon,colat,P_lm_spa2sph,w);
            sdeli_00 = deli_lm(1) - deli_00_prev;
            
            if include_lakes == 'y'
                % determine the depression adjacent to ice sheets;
                delP_j = calc_lake(ice_j_corr,oc_j,topo_j,lat_out,lon_out);
                delP_lm = spa2sph(delP_j,maxdeg,lon,colat,P_lm_spa2sph,w);
            else
                delP_lm = zeros(size(deli_lm));
            end

            % initial values for convergence
            conv = 'not converged yet';

            for k = 1:k_max % loop for sea level and topography iteration

                switch conv

                    case 'converged!'

                    case 'not converged yet'
                    
                    % calculate topography correction
                    TO = topo_0.*(oc_j-oc_0);
                    % expand TO function into spherical harmonics
                    TO_lm = spa2sph(TO,maxdeg,lon,colat,P_lm_spa2sph,w);
                    
                    % set up initial guess for sea level change
                    if k == 1 && topo_it == 1  % && t_it == 1
                        % initial guess of sea level change is just to distribute the
                        % ice over the oceans
                        sdelS_lm(t_it,:) = ocj_lm/ocj_lm(1)*(-rho_ice/rho_water*sdeli_00 + ...
                            (TO_lm(1) - TO_lm_prev(1)) - delP_lm(1)) - ...
                            (TO_lm - TO_lm_prev);
                    %elseif k == 1 && topo_it == 1 
                    %    sdelS_lm(t_it,:) = sdelS_lm(t_it-1,:);    
                    end
                    
                    delS_lm = delS_lm_prev + sdelS_lm(t_it,:);

                    % calculate change in loading
                    % delL is total change in loading
                    delL_lm = rho_ice*deli_lm + rho_water*delS_lm ...
                        + rho_sed*delSed_lm + rho_water*delP_lm;
                    % sdelL (small delta L) is incremental change in load -
                    % relative to last time step
                    sdelL_lm(t_it-1,:) = delL_lm - delL_lm_prev;
                    
                    % include rotation
                    % calculate degree two tidal k love number
                    
%                     % calculate contribution from rotation
%                     delLa_lm = calc_rot(delL_lm{t_it},k_L,k_T);
%                     sdelLa_lm = delLa_lm - delLa_lm_prev;


                    % calculate viscous contribution
                    
                    % beta contains the viscous love numbers for time t_it,
                    % row index goes over the time increments, column
                    % index goes over lm
                    if t_it == 2
                        V_lm = zeros(size(T_lm));
                    else
                        visc = beta_lm{t_it-1} .* sdelL_lm(1:t_it-2,:); 
                        % sum over all time increments. 
                        V_lm = sum(visc,1);
                    end
                    
                    
%                     % viscous rotation
%                     visc_T = squeeze(beta_T_lm(t_it-1,:,:)) .* sdelLa_lm; 
%                     % sum over all time increments.
%                     V_lm_T = sum(visc_T,1);


                    % calculate sea level perturbation
                    % add ice and sea level and multiply with love numbers
                    % DT doesn't load!
                    if include_rotation == 'y'
                        delSLcurl_lm_fl = E_lm .* T_lm .* delL_lm + T_lm .* V_lm + ...
                            1/g*E_lm_T.*delLa_lm + 1/g*V_lm_T + ...
                            T_lm .* V_lm;
                    else
                        delSLcurl_lm_fl = E_lm .* T_lm .* delL_lm + ...
                            T_lm .* V_lm;
                    end

                    % convert to spherical harmonics and subtract terms that are part
                    % of the topography to get the 'pure' sea level change
                    delSLcurl_fl = sph2spa(delSLcurl_lm_fl,maxdeg,lon,colat,P_lm_sph2spa);
                    delSLcurl = delSLcurl_fl - del_ice_corrected - del_DT - del_sed;


                    % compute and decompose RO
                    RO = delSLcurl.*oc_j;
                    RO_lm = spa2sph(RO,maxdeg,lon,colat,P_lm_spa2sph,w);

                    % calculate eustatic sea level perturbation (delta Phi / g)
                    delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*deli_lm(1) ...
                        - RO_lm(1) + TO_lm(1) - delP_lm(1));

                    
                    sdelS_lm_new = RO_lm + delPhi_g.*ocj_lm - TO_lm ...
                        - delS_lm_prev;
                    
                    
                    % calculate convergence criterion chi
                    chi = abs( (sum(abs(sdelS_lm_new)) - sum(abs(sdelS_lm(t_it,:)))) / ...
                        sum(abs(sdelS_lm(t_it,:))) );


                    % check convergence against the value epsilon
                    % If converged, set the variable conv to 'converged!' so that the
                    % calculation exits the loop. If not converged iterate again.
                    if chi < epsilon;
                        conv = 'converged!';
                        disp(['Finished time ' num2str(ice_time_new(t_it))...
                        'kyr. Number of iterations ' num2str(k) '. delphi is ' num2str(delPhi_g) ...
                        '. Lakes vol is ' num2str(delP_lm(1))])
                       % disp(['Converged after iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    else
                        conv = 'not converged yet';
                        %disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    end

                    % update sea sea surface height
                    sdelS_lm(t_it,:) = sdelS_lm_new;
               
                end

            end
            
            delS_lm_prev = delS_lm;
            TO_lm_prev = TO_lm;
            delL_lm_prev = delL_lm;
%            delLa_lm_prev = delLa_lm;
            deli_00_prev = deli_lm(1);
            
            % calculate overall perturbation of sea level over oceans
            % (spatially varying field and constant offset)
            delSL = delSLcurl + delPhi_g;
            
            % write in topography for next iteration
            topo(:,:,t_it) = - delSL + topo_0;
            
            if include_lakes == 'y'
                Lakes(:,:,t_it) = P_j;
            end

        end
        
        
        topo_diff = max(max(abs(topo(:,:,end) - topo0)));
        
        if topo_diff < max_topo_diff;
            conv_topo = 'converged!';
            disp(['Converged!! Number of topo iterations ' num2str(topo_it) ...
                '. Topo_diff is ' num2str(topo_diff)])
        else
            conv_topo = 'not converged yet';
            disp(['Not converged. Number of topo iterations ' num2str(topo_it) ...
                '. Topo_diff is ' num2str(topo_diff)])
        end

    end
    
    % update initial topography
    topo_initial(:,:,topo_it+1) = topo0 - (topo(:,:,end) - topo(:,:,1));
    
end

% calculate relative sea level (note that topography includes ice and we
% therefore need to subtract it here)
RSL = zeros(size(topo));
for i = 1:length(ice_time_new)
    RSL(:,:,i) = (topo(:,:,1) - ice_corrected(:,:,1)) - ...
        (topo(:,:,i) - ice_corrected(:,:,i));
end
%toc

topo_out = zeros(size(topo));
for i = 1:length(ice_time_new)
    topo_out(:,:,i) = topo0 - RSL(:,:,i);
end
%toc

%%

%file_path = '/Users/jackyaustermann/Dropbox/work_transfer/LakeBoneville/';
%lake_outline = importdata([file_path 'BonnevilleOutline.csv']);

figure
hold on

pcolor(lon_out,lat_out,RSL(:,:,end))
shading flat
axis([240 255 33 47])
colorbar

plot(lake_outline.data(:,1)+360,lake_outline.data(:,2),'ko')

% for i = 1:topo_it
%     topo_it_lm = spa2sph(topo_initial(:,:,i),maxdeg,lon,colat,P_lm_spa2sph,w);
%     topo_val(i) = sum(abs(topo_it_lm));
% end
% 
% chi_1 = abs((topo_val(2) - topo_val(1)) / topo_val(1));
% chi_2 = abs((topo_val(3) - topo_val(2)) / topo_val(2));

%% Plot results

% We only want the sea level change cause by melted ice, so subtract
% del_ice
fig_time = 21;
ind = find(ice_time_new==fig_time);
oc_1 = sign_01(topo(:,:,ind));
oc_2 = sign_01(topo(:,:,ind+1));

%plotSL = squeeze(topo_initial(:,:,2) - topo_initial(:,:,3));
plotSL = squeeze(RSL(:,:,ind) - RSL_b(:,:,ind));
% plotSL = squeeze(topo(:,:,2) - topo(:,:,1));
% plotSL = squeeze(ice_corrected(:,:,2) - ice_corrected(:,:,1));
%plotSL = squeeze(ice_area);
%plotSL = topo{ice_time_new==fig_time} - ice{ice_time_new==fig_time};
%plotSL = ice{ice_time_new==fig_time} - ice{ice_time_new==fig_time+1};

% plot
figure
m_proj('robinson','clongitude',0);
m_pcolor([lon_out(:,end/2+1:end)-360 lon_out(:,1:end/2)],lat_out,...
    [plotSL(:,end/2+1:end) plotSL(:,1:end/2)])
m_coast('color',[0 0 0]);
m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
colorbar
colormap(jet)

%% compare to Sam's results

% load bench_results/ICE5G_VM2_noTPW_floating.mat
% load bench_results/ICE5G_VM2_noTPW_nofloating.mat
% load ice5g_gl
plot_time = 21;

plotSL = squeeze(sl(:,:,times == plot_time));

% N = 256;
% [x,w] = GaussQuad(N);
% x_GL = acos(x)*180/pi - 90;
% lon_GL = linspace(0,360,2*N+1);
% lon_GL = lon_GL(1:end-1);
% [lon_tp,lat_tp] = meshgrid(lon_GL,x_GL);

% plot
figure
m_proj('robinson','clongitude',0);
m_pcolor([lon_out(:,end/2+1:end)-360 lon_out(:,1:end/2)],lat_out,...
    [plotSL(:,end/2+1:end) plotSL(:,1:end/2)])
m_coast('color',[0 0 0]);
m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
colorbar
colormap(jet)


%% Plot difference

% We only want the sea level change cause by melted ice, so subtract
% del_ice
fig_time = 21;
plotSL_JA = squeeze(RSL(:,:,ice_time_new==fig_time));
plotSL_SG = squeeze(sl(:,:,ice_time_new == fig_time));

plotSL = plotSL_JA - plotSL_SG;

% plot 
figure
m_proj('robinson','clongitude',0);
m_pcolor([lon_out(:,end/2+1:end)-360 lon_out(:,1:end/2)],lat_out,...
    [plotSL(:,end/2+1:end) plotSL(:,1:end/2)])
m_coast('color',[0 0 0]);
m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
colorbar
colormap(jet)

% caxis([-100 100])
%% Plot one time slice with ice and lakes

if include_lakes == 'y'

addpath ~/Documents/MATLAB/othercolor

time_fig = 12;

out = cptcmap('GMT_globe');

% We only want the sea level change cause by melted ice, so subtract
% del_ice
plotSL = topo{ice_time_new==time_fig} - ice{ice_time_new==time_fig};

% plot
figure
ax1 = axes;
%m_proj('lambert','longitude',[200 330],'latitude',[15 80]);
pcolor(ax1,lon_out,lat_out, plotSL)
%m_coast('color',[0 0 0]);
%m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
caxis([-6000 6000])
%colorbar
hold on

contour(ax1,lon_out,lat_out,plotSL,[0 0],'k')

% add ice plot
ice_plot = ice{ice_time_new==time_fig};
ice_plot(ice_plot == 0) = NaN;


ax2 = axes;
pcolor(ax2,lon_out,lat_out,ice_plot)
shading flat
caxis([-2000,5000])

% add lake plot plot
lake_plot = Lakes{ice_time_new==time_fig};
lake_plot(lake_plot == 0) = NaN;

ax3 = axes;
pcolor(ax3,lon_out,lat_out,lake_plot)
shading flat
caxis([0 700])

linkaxes([ax1,ax2,ax3])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

colormap(ax1,out)
colormap(ax2,flipud(othercolor('Blues4')))
colormap(ax3,othercolor('Purples9'))

axis([200 330 15 80])

end

%% Make video

make_vid = 1;

if make_vid == 1
else

% addpath ~/Documents/MATLAB/othercolor


clear A
out_g = cptcmap('GMT_globe');
out_r = cptcmap('GMT_relief');
out = [out_r(1:end/2,:); out_g(end/2+1:end,:)];

%axis([-0.6 0.6 -0.4 0.45])                    
%set(gca,'nextplot','replacechildren');

video_ind = 1;
init_ind = 1;

for i = find(ice_time_new == 22):find(ice_time_new == 0)
    
    if init_ind == 1
        A(1:length(find(ice_time_new == 1):length(ice_time_new))) ...
            = struct('cdata', [],'colormap', []);
        init_ind = 2;
    end
    
    figure(1);
    %hold on
    set(gcf, 'Position', [100 100 900 600])
    set(gcf,'color','w');
   

% We only want the sea level change cause by melted ice, so subtract
% del_ice
plotSL = topo{i} - ice{i};

% plot
ax1 = axes;
%m_proj('lambert','longitude',[200 330],'latitude',[15 80]);
pcolor(ax1,lon_out,lat_out, plotSL)
%m_coast('color',[0 0 0]);
%m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
caxis([-6000 6000])
%colorbar
hold on

contour(ax1,lon_out,lat_out,plotSL,[0 0],'k')

% add ice plot
ice_plot = ice{i};
ice_plot(ice_plot == 0) = NaN;

ax2 = axes;
pcolor(ax2,lon_out,lat_out,ice_plot/1000)
shading flat
caxis([-2,5])

% add lake plot plot
lake_plot = Lakes{i};
lake_plot(lake_plot == 0) = NaN;

ax3 = axes;
pcolor(ax3,lon_out,lat_out,lake_plot)
shading flat
caxis([0,500])

linkaxes([ax1,ax2,ax3])

ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

colormap(ax1,out)
colormap(ax2,flipud(othercolor('Blues4')))
colormap(ax3,othercolor('Purples9'))

axis([200 330 15 80])

text(215,25,[num2str(ice_time_new(i)) ' ka'],'FontSize',27,'color',[1 1 1])

h1 = colorbar(ax3,'Location','manual','position',[0.08 0.1 0.03 0.83]);
ylabel(h1,'lake depth (m)','FontSize',20)

% h2 = colorbar(ax3,'Location','manual','position',[0.91 0.1 0.03 0.83]);
% ylabel(h2,'lake depth (m)','FontSize',20)

    
%     title(['Time: ' num2str(ice_time_new(i)) ' ka cal BP'],...
%     'FontName','Arial','FontSize',12);
    
    pause(1.5);
    A(:,video_ind)=getframe(gcf); 
    video_ind = video_ind+1;
    close(1)
    
end

cd output_movies

movie2avi(A,'Lakes_video_512.avi', 'compression', 'None','fps',1,'quality',100)

cd ..

end