% Code to solve the elastic sea level equation following 
% Kendall et al., 2005 and Austermann et al., 2015

% J. Austermann 2015

% add paths when run for the first time.
% addpath SLFunctions
% addpath SLFunctions/mtimesx
% addpath '/Users/jackyaustermann/Documents/MATLAB/m_map'

%% Parameters & Input 
% Specify maximum degree to which spherical transformations should be done
maxdeg = 128;

% Some options to choose from
include_lakes = 'n'; % choose between y (for yes) and n (for no)
include_rotation = 'n'; % choose between y (for yes) and n (for no)
include_ice_check = 'y'; % choose between y (for yes) and n (for no)
include_sediments = 'n'; % choose between y (for yes) and n (for no)
include_dynamic_topography = 'n'; % choose between y (for yes) and n (for no)

% parameters
rho_ice = 920;
rho_water = 1000;
rho_sed = 2300;
g = 9.80665;


% The following steps help speed up the calculations
% Set up Gauss Legendre grid onto which to interpolate all grids
N = maxdeg; %change size if you want to run it at higher resolution than the maxdeg
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

ice = single(zeros(length(x_GL),length(lon_GL),length(ice_time)));
ice_time_new = zeros(size(ice_time));

ice_lat = [90; ice_lat; -90];
ice_long = [ice_long, 360];

for i = 1:length(ice_time)

    % put onto on Gauss Legendre grid
    ice_nointerp = squeeze(ice5g_grid(i,:,:));
    
    % add rows at top and bottom, and right
    ice_extended = [zeros(1,length(ice_long)-1); ice_nointerp; ...
        ice_nointerp(1,end)*ones(1,length(ice_long)-1)];
    
    ice_extended_2 = [ice_extended, ice_extended(:,1)];

    % interpolate ice masks on Gauss Legendre grid
    ice_interp = interp2(ice_long,ice_lat,ice_extended_2,lon_out, lat_out);
    
    ice(:,:,length(ice_time)-i+1) = ice_interp;
    ice_time_new(i) = ice_time(length(ice_time)-i+1);
end


% --------------------------------
% TOPOGRAPHY
% --------------------------------

% load preloaded etopo, which includes interpolated fields onto different
% sized Gauss Legendre Grids (hence avoiding interpolating twice)
load topo_SL

% interpolate topography grid onto Gauss Legendre Grid
if N == 64
    topo_pres = topo_bed_64 + ice(:,:,end);
elseif N == 128
    topo_pres = topo_bed_128 + ice(:,:,end);
elseif N == 256
    topo_pres = topo_bed_256 + ice(:,:,end);
elseif N == 512
    topo_pres = topo_bed_512 + ice(:,:,end);
elseif N == 1024
    topo_pres = topo_bed_1024 + ice(:,:,end);
else
    topo_pres = interp2(lon_topo,lat_topo,topo_bed,lon_out,lat_out) + ice(:,:,end);
end


% --------------------------------
% DYNAMIC TOPOGRAPHY
% --------------------------------

%del_DT = zeros(size(lon_out));
DT = single(zeros(length(x_GL),length(lon_GL),length(ice_time_new)));

% --------------------------------
% SEDIMENT
% --------------------------------

% sediment history 
sed = single(zeros(length(x_GL),length(lon_GL),length(ice_time_new)));



%% Set up love number input

% prepare love numbers in suitable format and calculate T_lm and E_lm 
% to calculate the fluid case, switch h_el to h_fl, k_el to k_fl and same
% for tidal love numbers
load SavedLN/prem.l90C.umVM2.lmVM2.mat
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
h_lm_tide = love_lm(h_el_tide,maxdeg);
k_lm_tide = love_lm(k_el_tide,maxdeg);

E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);

E_lm_T = 1 + k_lm_tide - h_lm_tide;


% calculate betas for viscous contribution
beta_l = cell(length(ice_time_new)-1);

for t_it = 2:length(ice_time_new)
    
    for n = 2:t_it-1
        
        beta = zeros(maxdeg, 1);
        for lm = 1:maxdeg
            num_mod = mode_found(lm);
            beta(lm) = sum((k_amp(lm,1:num_mod) - h_amp(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n)))));
        end
        
        beta_l{t_it-1}(n-1,:) = [0; beta]; % add 0,0 LN

    end
end

% Since the Love Numbers only depend on l, but not on m we need to compute
% a counter that maps from a _l vector to a _lm vector
beta_counter = ones(size(h_lm));
l_it = 1;
for lm_it = 1:length(h_lm)
    if lm_it == l_it*(l_it+1)/2
        beta_counter(lm_it+1) = beta_counter(lm_it)+1;
        l_it = l_it+1;
    else
        beta_counter(lm_it+1) = beta_counter(lm_it);
    end
end


%% Solve sea level equation (after Kendall 2005, Dalca 2013 & Austermann et al. 2015)
tic
k_max = 10;   % maximum number of iterations
epsilon = 10^-4; % convergence criterion

topo_it_max = 1;   % maximum number of iterations
max_topo_diff = 10; % convergence criterion

% 0 = before
% j = after

% set up initial topography
topo_initial = zeros(length(x_GL),length(lon_GL),topo_it_max+1);
topo_initial(:,:,1) = topo_pres - ice(:,:,end) + ice(:,:,1); % already includes ice, DT and sediments

% initial topography with the guess that topography is the same as present at every
% point in time; topography is a 3D array; access topography at time x
% like this topo(:,:,x) [or for plotting squeeze(topo(:,:,x))]
topo = zeros(length(x_GL),length(lon_GL),length(ice_time_new));
for i = 2:length(ice_time_new)
    topo(:,:,i) = topo_pres - ice(:,:,end) + ice(:,:,i);
end

% initialize sediment and DT changes
delSed_lm = zeros(length(ice_time_new)-1,length(h_lm));
for t_it = 2:length(ice_time_new) 
    % calculate change in sediments and decompose into spherical harmonics
    del_sed = sed(:,:,t_it) - sed(:,:,1);
    if include_sediments == 'y'
        delSed_lm(t_it-1,:) = spa2sph(del_sed,maxdeg,lon,colat,P_lm_spa2sph,w);
    else
        delSed_lm(t_it-1,:) = zeros(size(h_lm));
    end
end


% initialize Lakes, ocean load and corrected ice
Lakes = zeros(size(topo));
sdelS_lm = zeros(length(ice_time_new),length(h_lm));
ice_corrected = ice;


% initial values for convergence
conv_topo = 'not converged yet';


% TOPOGRAPHY ITERATION
for topo_it = 1:topo_it_max;
    
    switch conv_topo

        case 'converged!'

        case 'not converged yet'
            
        % initialize arrays for each timestep
        delL_lm_prev = zeros(1,length(h_lm));
        delS_lm_prev = zeros(1,length(h_lm));
        TO_lm_prev = zeros(1,length(h_lm));
        % delLa_lm_prev = zeros(1,length(h_lm)); for rotation
        deli_00_prev = 0;
        delP_00_prev = 0;
        
        sdelL_lm = zeros(length(ice_time_new)-1,length(h_lm));

        % update new initial topography
        topo(:,:,1) = topo_initial(:,:,topo_it);

        % remove the corrected ice model and add initial ice model back on
        % this needs to be done to calculate the updated corrected ice model
        % (only necessary for later iterations, not the first one)
        if topo_it == 1
        else
            for i = 1:length(ice_time_new) 
                topo(:,:,i) = topo(:,:,i) - ice_corrected(:,:,i) + ice(:,:,i);
            end
        end

        % recompute corrected ice model
        % do grounded ice check to calculate the corrected ice model
        for i = 1:length(ice_time_new)
            if include_ice_check == 'y'
                 % check ice model for floating ice
                 check1 = sign_01(-topo(:,:,i) + ice(:,:,i));
                 check2 = sign_01(+topo(:,:,i) - ice(:,:,i)) .* ...
                     (sign_01(-ice(:,:,i)*rho_ice - (topo(:,:,i) - ice(:,:,i))*rho_water));

                 ice_corrected(:,:,i) = check1.*ice(:,:,i) + check2.*ice(:,:,i);
            else
                % if the floating ice check is set to 'n', don't change the
                % ice model
                 ice_corrected(:,:,i) = ice_j;
            end
        end

        % update all topographies with the new / corrected ice model
        for i = 1:length(ice_time_new) 
            topo(:,:,i) = topo(:,:,i) - ice(:,:,i) + ice_corrected(:,:,i);
        end

        % assign topography of time 0 and calculate the ocean function
        topo_0 = topo(:,:,1);
        oc_0 = sign_01(topo_0);
        oc0_lm = spa2sph(oc_0,maxdeg,lon,colat,P_lm_spa2sph,w);
        ocj_lm_prev = oc0_lm;
        

        % TIME ITERATION
        for t_it = 2:length(ice_time_new) 

            % Assign topography and ocean function of time t_it to the
            % index j
            topo_j = topo(:,:,t_it);
            oc_j = sign_01(topo_j);
            ocj_lm = spa2sph(oc_j,maxdeg,lon,colat,P_lm_spa2sph,w);

            % calculate topography correction
            TO = topo_0.*(oc_j-oc_0);
            TO_lm = spa2sph(TO,maxdeg,lon,colat,P_lm_spa2sph,w);
            
            % calculate change in sediments and dynamic topography
            del_sed = sed(:,:,t_it) - sed(:,:,1);
            del_DT = DT(:,:,t_it) - DT(:,:,1);
            
            % calculate the change in ice model
            del_ice_corrected = ice_corrected(:,:,t_it) - ice_corrected(:,:,1);
            deli_lm = spa2sph(del_ice_corrected,maxdeg,lon,colat,P_lm_spa2sph,w);
            % calculate the incremental increase in ice volume
            sdeli_00 = deli_lm(1) - deli_00_prev;


            % When including lakes, calculate their location and size
            if include_lakes == 'y'
                % determine the depression adjacent to ice sheets;
                delP_j = calc_lake(ice_j_corr,oc_j,topo_j,lat_out,lon_out);
                delP_lm = spa2sph(delP_j,maxdeg,lon,colat,P_lm_spa2sph,w);
            else
                delP_lm = zeros(size(deli_lm));
            end


            % initial values for convergence
            conv = 'not converged yet';

            % SEA LEVEL EQUATION ITERATION
            for k = 1:k_max 

                switch conv

                    case 'converged!'

                    case 'not converged yet'

                    % set up initial guess for sea level change
                    if k == 1 && topo_it == 1
                        % initial guess of sea level change is just to distribute the
                        % ice over the oceans
                        % use slightly different initial guess than Kendall

                        sdelS_lm(t_it,:) = ocj_lm_prev/ocj_lm_prev(1)*...
                            (-rho_ice/rho_water*sdeli_00 + ...
                            TO_lm(1)-TO_lm_prev(1) - (delP_lm(1)-delP_00_prev)) ...
                            - TO_lm - TO_lm_prev;
                    end

                    % calculate total change in ocean load
                    delS_lm = delS_lm_prev + sdelS_lm(t_it,:);

                    % calculate change in loading
                    % delL is total change in loading
                    delL_lm = rho_ice*deli_lm + rho_water*delS_lm ...
                        + rho_sed*delSed_lm(t_it-1,:) + rho_water*delP_lm;
                    % sdelL (small delta L) is incremental change in load -
                    % relative to last time step
                    sdelL_lm(t_it-1,:) = delL_lm - delL_lm_prev;


                    % calculate viscous contribution to deformation
                    % beta contains the viscous love numbers for time t_it,
                    % row index goes over the time increments, column
                    % index goes over lm
                    if t_it == 2
                        V_lm = zeros(size(T_lm));
                    else
                        for lm_it = 1:length(h_lm)
                            V_lm(lm_it) = beta_l{t_it-1}(:,beta_counter(lm_it))'...
                                * sdelL_lm(1:t_it-2,lm_it);
                        end
                    end

                    
                    % include rotation
                    % calculate degree two tidal k love number

                    % TBD

                    
                    % calculate sea level perturbation
                    if include_rotation == 'y'
                        disp('This code doesnt do rotation yet')
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


                    % calculate incremental change in sea surface height
                    sdelS_lm_new = RO_lm + delPhi_g.*ocj_lm - TO_lm ...
                        - delS_lm_prev;


                    % calculate convergence criterion chi
                    chi = abs((sum(abs(sdelS_lm_new)) - sum(abs(sdelS_lm(t_it,:)))) / ...
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
                    elseif chi < epsilon && k == k_max;
                        conv = 'not converged yet';
                        disp(['Finished time ' num2str(ice_time_new(t_it))...
                        'kyr. Run has not converged. Chi is  ' num2str(chi)])
                    else
                        conv = 'not converged yet';
                        %disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    end

                    % update sea sea surface height
                    sdelS_lm(t_it,:) = sdelS_lm_new;

                end

            end

            % update values
            delS_lm_prev = delS_lm;
            TO_lm_prev = TO_lm;
            delL_lm_prev = delL_lm;
    %            delLa_lm_prev = delLa_lm; % for rotation
            deli_00_prev = deli_lm(1);
            delP_00_prev = delP_lm(1);

            % calculate overall perturbation of sea level over oceans
            % (spatially varying field and constant offset)
            delSL = delSLcurl + delPhi_g;

            % calculate topography and use as initial guess for next
            % iteration
            topo(:,:,t_it) = - delSL + topo_0;
            ocj_lm_prev = ocj_lm;

            if include_lakes == 'y'
                Lakes(:,:,t_it) = P_j;
            end

        end

        % correct the present-day topography by the corrected ice model, so
        % that we're comparing two topographies with the same ice model
        topo_pres_ice_corrected = topo_pres - ice(:,:,end) + ice_corrected(:,:,end);

        % calculate the maximum difference between the two
        topo_diff = max(max(abs(topo(:,:,end) - topo_pres_ice_corrected)));

        % check convergence
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
    topo_initial(:,:,topo_it+1) = topo_pres_ice_corrected - (topo(:,:,end) - topo(:,:,1));
    
end

% calculate relative sea level (note that topography includes ice and we
% therefore need to subtract it here)
RSL = zeros(size(topo));
for i = 1:length(ice_time_new)
    RSL(:,:,i) = (topo(:,:,end) - ice_corrected(:,:,end)) - ...
        (topo(:,:,i) - ice_corrected(:,:,i));
end
toc


%% Plot results

% We only want the sea level change cause by melted ice, so subtract
% del_ice
fig_time = 21;
ind = find(ice_time_new==fig_time);

plotSL = squeeze(RSL(:,:,ind) - RSL_old(:,:,ind));

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