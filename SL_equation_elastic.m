% Code to solve the elastic sea level equation following 
% Kendall et al., 2005 and Austermann et al., 2015

% J. Austermann 2015

% add paths when run for the first time.
% addpath SLFunctions
% addpath '/Users/jackyaustermann/Documents/MATLAB/m_map'

%% Parameters & Input 
% Specify maximum degree to which spherical transformations should be done
maxdeg = 512;

% parameters
rho_ice = 916;
rho_water = 1000;
rho_sed = 2300;
g = 9.81;


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


% --------------------------------
% ICE
% --------------------------------

% load WAIS 
load ice5g_griddata
ice = cell(size(ice_time));
ice_time_new = zeros(size(ice_time));

for i = 1:length(ice_time)

    % already on Gauss Legendre grid
    ice_nointerp = squeeze(ice5g_grid(i,:,:));

    % interpolate ice masks on common grid
    ice_interp = interp2(ice_long,ice_lat,ice_nointerp,lon_out, lat_out);

    % patch in zeros
    ice_interp(isnan(ice_interp) == 1) = 0;
    
    ice{length(ice_time)-i+1} = ice_interp;
    ice_time_new(i) = ice_time(length(ice_time)-i+1);
end

%del_ice = ice_j - ice_0; 


% --------------------------------
% DYNAMIC TOPOGRAPHY
% --------------------------------

del_DT = zeros(size(lon_out));


% --------------------------------
% SEDIMENT
% --------------------------------

del_sed = zeros(size(lon_out));


% --------------------------------
% TOPOGRAPHY
% --------------------------------

% load preloaded etopo (including ice) as topo_orig, lon_topo, lat_topo
load topo_SL

% interpolate topography grid onto Gauss Legendre Grid
topo0 = interp2(lon_topo,lat_topo,topo_orig,lon_out, lat_out);



%% Set up love number input

% prepare love numbers in suitable format and calculate T_lm and E_lm 
% to calculate the fluid case, switch h_el to h_fl, k_el to k_fl and same
% for tidal love numbers
load SavedLN/LN_l120_ump5_lm5_512
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
h_lm_tide = love_lm(h_el_tide,maxdeg);
k_lm_tide = love_lm(k_el_tide,maxdeg);

E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);

E_lm_T = 1 + k_lm_tide - h_lm_tide;

% can switch this in if you want to exclude rotational effects
% E_lm_T = zeros(size(E_lm_T));

%% Solve sea level equation (after Kendall 2005, Dalca 2013 & Austermann et al. 2015)
tic
k_max = 10;   % maximum number of iterations
epsilon = 10^-3; % convergence criterion

topo_it_max = 5;   % maximum number of iterations
max_topo_diff = 10; % convergence criterion

% 0 = before
% j = after

% set up initial topography and ocean function
topo_initial{1} = topo0 + ice{1}; % already includes ice and dynamic topography

% no proglacial lakes in step 0 
P_0 = zeros(size(topo0));
Lakes{1} = zeros(size(topo0));

topo = cell(1,length(ice_time));
SL_change = cell(1,length(ice_time));
Lakes = cell(1,length(ice_time));
deli_lm = cell(1,length(ice_time));
delS_lm = cell(1,length(ice_time));
delP_lm = cell(1,length(ice_time));

deli_lm{1} = zeros(size(h_lm));
delS_lm{1} = zeros(size(h_lm));
delP_lm{1} = zeros(size(h_lm));

V_lm_0 = zeros(size(h_lm));


% initial values for convergence
conv_topo = 'not converged yet';

for topo_it = 1:topo_it_max;
    
    topo_0 = topo_initial{topo_it};
    oc_0 = sign_01(topo_0);

    % expand ocean function into spherical harmonics
    oc0_lm = spa2sph(oc_0,maxdeg,lon,colat,P_lm_spa2sph);
    
    switch conv_topo

        case 'converged!'

        case 'not converged yet'

        for t_it = 2:length(ice_time) % loop over time

            visc = zeros(t_it-1,length(h_lm));

            topo{t_it-1} = topo_0;

            ice_0 = ice{t_it-1};
            ice_j = ice{t_it};

            % set up initial guess for new topography
            topo_j = topo{t_it-1} + ice_j-ice_0; % del_ice is negative -> subtract ice that is melted
            oc_j = sign_01(topo_j);

            % calculate change in sediments and decompose into spherical harmonics
            % set to zero as is but can be set as the change between the sed load. 
            Sed_lm = zeros(size(oc0_lm));%spa2sph(del_sed,maxdeg,lon,colat,P_lm_spa2sph);

            % initial values for convergence
            conv = 'not converged yet';

            for k = 1:k_max % loop for sea level and topography iteration

                switch conv

                    case 'converged!'

                    case 'not converged yet'

                    % expand ocean function into spherical harmonics
                    ocj_lm = spa2sph(oc_j,maxdeg,lon,colat,P_lm_spa2sph);


                    % CHECK ICE MODEL 
                    % check ice model for floating ice
                    % check1 = sign_01(-topo_j + ice_j);
                    % check2 = sign_01(+topo_j - ice_j) .* ...
                    %    (sign_01(-ice_j*rho_ice - (topo_j - ice_j)*rho_water));

                    % ice_j_corr = check1.*ice_j + check2.*ice_j;
                    ice_j_corr = ice_j;
                    del_ice_corrected = ice_j_corr - ice_0; 

                    deli_lm{t_it} = spa2sph(del_ice_corrected,maxdeg,lon,colat,P_lm_spa2sph);

                    % determine the depression adjacent to ice sheets;
                    P_j = calc_lake(ice_j_corr,oc_j,topo_j,lat_out,lon_out);
                    delP = P_j - P_0;
                    delP_lm{t_it} = spa2sph(delP,maxdeg,lon,colat,P_lm_spa2sph);

                    % calculate topography correction
                    TO = topo_0.*(oc_j-oc_0);
                    % expand TO function into spherical harmonics
                    TO_lm = spa2sph(TO,maxdeg,lon,colat,P_lm_spa2sph);


                    % set up initial guess for sea level change
                    if k == 1
                        % initial guess of sea level change is just to distribute the
                        % ice over the oceans
                        delS_lm{t_it} = ocj_lm/ocj_lm(1)*(-rho_ice/rho_water*deli_lm{t_it}(1) + ...
                            TO_lm(1) - delP_lm{t_it}(1));
                        % convert into spherical harmonics
                        delS_init = sph2spa(delS_lm{t_it},maxdeg,lon,colat,P_lm_sph2spa);

                    end

                    % calculate elastic loading term
                    L_lm = rho_ice*deli_lm{t_it} + rho_water*delS_lm{t_it} + rho_sed*Sed_lm + ...
                        rho_water*delP_lm{t_it};

                    % calculate contribution from rotation
                    La_lm = calc_rot(L_lm,k_el,k_el_tide);


                    % calculate viscous contribution
                    for n = 1:t_it-1
                        beta = zeros(maxdeg, 1);
                        for lm = 1:maxdeg
                            num_mod = mode_found(lm);
                            beta(lm) = sum((k_amp(lm,1:num_mod) - h_amp(lm,1:num_mod)) ...
                                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                                * (ice_time(t_it) - ice_time(n)))));
                        end
                        beta_lm = love_lm(beta, maxdeg);
                        visc(n, 1:length(beta_lm)) = beta_lm .* (rho_ice * deli_lm{n} ...
                            + rho_water*delS_lm{n} + rho_sed*Sed_lm + rho_water*delP_lm{n}); 
                    end

                    V_lm_j = sum(visc,1);
                    V_lm = V_lm_j - V_lm_0;

                    % calculate sea level perturbation
                    % add ice and sea level and multiply with love numbers
                    % DT doesn't load!
                    delSLcurl_lm_fl = E_lm .* T_lm .* L_lm + T_lm .* V_lm + ...
                        1/g*E_lm_T.*La_lm;
                    % THIS DOESNT YET INCLUDE VISCOUS ROTATION

                    % convert to spherical harmonics and subtract terms that are part
                    % of the topography to get the 'pure' sea level change
                    delSLcurl_fl = sph2spa(delSLcurl_lm_fl,maxdeg,lon,colat,P_lm_sph2spa);
                    delSLcurl = delSLcurl_fl - del_ice_corrected - del_DT - del_sed;


                    % compute and decompose RO
                    RO = delSLcurl.*oc_j;
                    RO_lm = spa2sph(RO,maxdeg,lon,colat,P_lm_spa2sph);

                    % calculate eustatic sea level perturbation (delta Phi / g)
                    delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*deli_lm{t_it}(1) ...
                        - RO_lm(1) + TO_lm(1) - delP_lm{t_it}(1));


                    % calculate overall perturbation of sea level over oceans
                    % (spatially varying field and constant offset)
                    delSL = delSLcurl + delPhi_g;


                    % update topography and ocean function
                    topo_j = - delSL + topo_0;
                    oc_j = sign_01(topo_j);


                    % calculate change in ocean height and decompose
                    delS_new = delSL.*oc_j -  topo_0.*(oc_j-oc_0);
                    delS_lm_new = spa2sph(delS_new,maxdeg,lon,colat,P_lm_spa2sph);


                    % calculate convergence criterion chi
                    chi = abs( (sum(abs(delS_lm_new)) - sum(abs(delS_lm{t_it}))) / ...
                        sum(abs(delS_lm{t_it})) );


                    % check convergence against the value epsilon
                    % If converged, set the variable conv to 'converged!' so that the
                    % calculation exits the loop. If not converged iterate again.
                    if chi < epsilon;
                        conv = 'converged!';
                        disp(['Finished time ' num2str(ice_time_new(t_it))...
                        'kyr. Number of iterations ' num2str(k) '. delphi is ' num2str(delPhi_g) ...
                        '. del Lakes vol is ' num2str(delP_lm{t_it}(1))])
                       % disp(['Converged after iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    else
                        conv = 'not converged yet';
                        %disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    end

                    % update sea sea surface height
                    delS_lm{t_it} = delS_lm_new;
                end

            end

            topo{t_it} = topo_j;
            topo_0 = topo_j; % initial topography for next timestep
            oc_0 = sign_01(topo_0);
            SL_change{t_it} = delSL + del_ice_corrected;
            P_0 = P_j;
            Lakes{t_it} = P_j;
            V_lm_0 = V_lm_j;

        end
        
        
        topo_diff = max(max(topo{end} - topo0));
        
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
    
    topo_initial{topo_it+1} = topo_initial{topo_it} + (topo0 - topo{end});
    
end
toc


%% Plot results

% We only want the sea level change cause by melted ice, so subtract
% del_ice
plotSL = SL_change{ice_time_new==40};

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

%% Make video


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