%
%  Copyright (C) 2016 - 2024 by J. Austermann.
%  Code is commented by J. Austermann & S. Chester
%  This file is part of SLcode.
%  SLcode is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%  SLcode is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  <http://www.gnu.org/licenses/>.


% Code to solve the elastic sea level equation following 
% Kendall et al., 2005 and Austermann et al., 2015

% Make sure these folders are located in your working directory (or add the
% complete path to files. 
addpath SLFunctions/
addpath SavedLN/

%% Parameters & Input 
% Set Spatial Resolution: Specify maximum degree to which spherical
% transformations should be calculated. (Degree 256 = ~80 km). Choose a
% value 2^n for computational efficiency, i.e. 64, 128, 256, 512
maxdeg = 256; 

% Include the effect of changes in Earth's Rotation:
include_rotation = 'y'; % y=yes and n=no

% Include a check that looks whether parts of ice sheets are floating or
% not. If they are, remove them from the load (since they are neutrally
% bouyant). 
include_ice_check = 'y'; % y=yes and n=no

% parameters
rho_ice = 920; % density of ice
rho_water = 1000; % density of water
g = 9.80665; % gravitational accelaration 

% Set up a Gauss Legendre (GL) grid onto which to interpolate all grids
% (Using a GL grid helps when transforming to spherical harmonics)
N = maxdeg; 
[x,w] = GaussQuad(N);
x_GL = acos(x)*180/pi - 90;
lon_GL = linspace(-180,180,2*N+1);
lon_GL = lon_GL(1:end-1);

[lon_out,lat_out] = meshgrid(lon_GL,x_GL); % x and y (lon and lat) matrices.

% Precompute legendre polynomials (used in spherical harmonic transform)
P_lm = cell(N+1,1); 
for l=0:N
    P_lm{l+1} = legendre(l,x,'norm');
end


%% -------------------------------
% ICE SHEET RECONSTRUCTION
% --------------------------------

% The ice history needs to have the form of a 3D matrix following the form
% [longitude, latitude, time]. The ice needs to be interpolated onto the
% Gauss Legendre grid constructed above. The name of the variable that
% stores the ice history should be 'ice'. 

% The time vector that goes with the above ice history should be a 1D
% vector called 'ice_time_new'. This determines the temporal resolution of
% the calculation. Time should be in ka, i.e. in thousands of years from
% past to present. For example ice_time_new = [21, 20, 19 ... 0] would be a
% run from the last glacial maximum (21 ka) to present (0 ka). 

%load fingerprint_melt.mat
%ice = WAIS_melt;
%ice = Greenland_melt;

load ice7g_122ka_GL256.mat
% - This is the ICE-7G model (Roy and Peltier, 2017; 2018). The published
% model extends from 26 ka to present. We extended this model to 122 ka using 
% the Waelbroeck et al. (2002) sea level reconstruction, assuming times of
% similar global mean sea level correspond to similar ice sheet geometry.
% - The model is already interpolated onto a degree 256 Gauss Legendre
% grid. 

% determine what index refers to present-day for the topography iteration
ind_pres = find(ice_time_new == 0);

%% --------------------------------
% TOPOGRAPHY
% --------------------------------

% load preloaded etopo, which includes interpolated fields onto different
% sized Gauss Legendre Grids (to avoid interpolating multiple times)
load topo_SL

% interpolate topography grid onto Gauss Legendre Grid
if N == 64
    topo_bed_64 = [topo_bed_64(:,end/2+1:end) topo_bed_64(:,1:end/2)];% input on a longitude grid from 0-360 but we want -180 to 180
    topo_pres = topo_bed_64 + ice(:,:,ind_pres);
elseif N == 128
    topo_bed_128 = [topo_bed_128(:,end/2+1:end) topo_bed_128(:,1:end/2)];
    topo_pres = topo_bed_128 + ice(:,:,ind_pres);
elseif N == 256
    topo_bed_256 = [topo_bed_256(:,end/2+1:end) topo_bed_256(:,1:end/2)];
    topo_pres = topo_bed_256 + ice(:,:,ind_pres);
elseif N == 512
    topo_bed_512 = [topo_bed_512(:,end/2+1:end) topo_bed_512(:,1:end/2)];
    topo_pres = topo_bed_512 + ice(:,:,ind_pres);
elseif N == 1024
    topo_bed_1024 = [topo_bed_1024(:,end/2+1:end) topo_bed_1024(:,1:end/2)];
    topo_pres = topo_bed_1024 + ice(:,:,ind_pres);
else
    topo_pres = interp2(lon_topo,lat_topo,topo_bed,lon_out,lat_out) + ice(:,:,ind_pres);
end
 
% calculate the 'ocean function' (1=ocean, 0=not ocean) for the present day. 
oc_pres = sign_01(topo_pres);
ocpres_lm = sphere_har(oc_pres,maxdeg,N,P_lm);
oc_area = ocpres_lm(1);


%% --------------------------------
% Set up love number input
% --------------------------------

% The Sea Level Equation is solved using Green's functions, which are
% described using Love Numbers. These love numbers represent the 
% response of the solid Earth and gravity field to surface perturbations 
% (ice load changes). See Kendall et al. (2005) & Peltier (1974) for details.

% '_tide' and '_T' indicates tidal love numbers which are used for calculating the
% rotation component. 

% Load love numbers in suitable format
% Love numbers are calculated using the collocation method (Mitrovica and
% Peltier, 1992). 
% 'prem' refers to the elastic/density structure (from Dziewonsku and Anderson, 1981)
% 'lXX' indicated a fully elastic lithosphere of XX thickness (km)
% 'umXX' indicates an upper mantle viscosity of XX *10^21 Pa s. 
% 'lmXX' indicates a lower mantle viscosity of XX *10^21 Pa s. 
% note that a 'p' indicates '0.' so 'p4' is '0.4'

load SavedLN/prem.l100C.ump4.lm4.mat %file containing the love numbers

% Read in love numbers and uce love_lm to order them appropriately. 
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
h_lm_tide = love_lm(h_el_tide,maxdeg);
k_lm_tide = love_lm(k_el_tide,maxdeg);

% The degree 1 component is in a reference frame that is shifting significantly. 
% We will change this here into a center of mass of the solid Earth reference 
% frame in which k_amp(1,:) is 0. 
h_amp(1,:) = h_amp(1,:) - k_amp(1,:);
k_amp(1,:) = zeros(size(k_amp(1,:)));

E_lm = 1 + k_lm - h_lm; %variable containting the total elastic deformation of the gravity field and solid Earth decribed by the love numbers. 
T_lm = get_tlm(maxdeg); %variable representing constants in the sea level equation (Earth's radius, mass, pi etc.)

E_lm_T = 1 + k_lm_tide - h_lm_tide; % components necessary for calculating changes in rotation

% Calculate 'betas'. These are functions used in the SLE to represent 
% the viscous deformation (See Kendall et al., 2005 for details)
beta_l = cell(length(ice_time_new)-1,1);
beta_konly_l = cell(length(ice_time_new)-1,1);

for t_it = 2:length(ice_time_new)
    for n = 2:t_it-1
        beta = zeros(maxdeg, 1);
        for lm = 1:maxdeg
            num_mod = mode_found(lm);
            beta(lm) = sum((k_amp(lm,1:num_mod) - h_amp(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n)))));
        end
        beta_l{t_it-1}(n-1,:) = [0; beta]; % add 0 LN
        % for rotation only needed for degree 2
        lm = 2;
        num_mod = mode_found(lm);
        beta_konly_l{t_it-1}(n-1) = sum((k_amp(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n)))));
    end
end


% calculate tidal betas (for rotation)
beta_tide = cell(length(ice_time_new)-1,1);
beta_konly_tide = cell(length(ice_time_new)-1,1);

for t_it = 2:length(ice_time_new)
    for n = 2:t_it-1
        beta = zeros(maxdeg, 1);
        for lm = 1:maxdeg
            num_mod = mode_found(lm);
            beta(lm) = sum((k_amp_tide(lm,1:num_mod) - h_amp_tide(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n))))); 
        end
        beta_tide{t_it-1}(n-1,:) = [0; beta]; % add 0 LN
        % for rotation only needed for degree 2
        lm = 2;
        num_mod = mode_found(lm);
        beta_konly_tide{t_it-1}(n-1) = sum((k_amp_tide(lm,1:num_mod)) ...
                ./spoles(lm,1:num_mod).* (1 - exp(- spoles(lm,1:num_mod) ...
                * (-ice_time_new(t_it) + ice_time_new(n)))));
    end
end

% Initiate mapping from l to lm. 
% The 'beta_counter' function lets us map the betas (2D cell variables) to
% the love number 1D vectors which are ordered by degree and order. (e.g.
% element 1 is degree0, order0, element 2 is degree1, order0, element 3 is
% degree1, order2, and so on). 
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


%% -------------------------------
% Solve sea level equation 
% --------------------------------
% (after Kendall 2005, Dalca et al. 2013 & Austermann et al. 2015)

% Solving the SLE involves three nested for loops. 
%  1) An outermost 'topography' iteration. Used to find the starting
%     topography that results in a prescribed final (present-day) topography.
%  2) A middle 'time interaiton'. Loops through time.
%  3) And an innermost 'sea level equation (SLE)' iteration, which iteratively
%     updates the ocean load. 

% Variable Notes: 
% 'ice_time_new' = time vector, ordered old to young.
% 'ice' = the ice model; 3D matrix: [lat, lon, time].
% 'ice_corrected' = the ice model after the marine based correction is made
% 'topo' = is a 3D vector; 3D matrix: [lat, lon, time].
% 'lon_out' = longitude matrix (corresponds to ice model, topography, and all code output)
% 'lat_out' = latitude matrix (same as lon_out)
% 'del_' = indiated cummulative changes/deformation 
% 'sdel_' = indicates incremental (single time-step) changes
% '_lm' = indicates that the variable is in the spherical harmonic form 
% 'SLcurl' = non-uniform changes in sea level (no uniform shift). 

% time syntax:
% 0 = before
% j = after

tic
k_max = 10;   % maximum number of 'SLE' iterations (innermost loop)
epsilon = 10^-4; % convergence criterion for the SLE loop

topo_it_max = 3;   % maximum number of 'Topography' iterations 
max_topo_diff = 1; % convergence criterion for topography loop

% set up initial (t=122ka) topography and ocean function
topo_initial = zeros(length(x_GL),length(lon_GL),topo_it_max+1);
% Topo_pres includes modern ice, so we subtract it here and add to the paleo ice sheet 
topo_initial(:,:,1) = topo_pres - ice(:,:,ind_pres) + ice(:,:,1); % 

% initial topography guess: topography is the same as present at every
% point in time; topography is a 3D vector; access topography at time x
% like this topo(:,:,x) [or for plotting squeeze(topo(:,:,x))]

% For the first time through the loop we guess that the initial topography 
% is the same as present at every point in time.
topo = zeros(length(x_GL),length(lon_GL),length(ice_time_new));
for i = 2:length(ice_time_new)
    topo(:,:,i) = topo_pres - ice(:,:,ind_pres) + ice(:,:,i);
end

% initialize more variables
ice_corrected = ice; 
sdelS_lm = zeros(length(ice_time_new),length(h_lm));

% initial values for convergence
conv_topo = 'not converged yet';

% TOPOGRAPHY ITERATION
for topo_it = 1:topo_it_max
    
    switch conv_topo

        case 'converged!' % end the loop!

        case 'not converged yet'
            
        % initialize for each time loop
        delL_lm_prev = zeros(1,length(h_lm)); %cummulative change in Load (ice and ocean mass)
        delS_lm_prev = zeros(1,length(h_lm)); %cummulative change in S (sea surface height) 
        TO_lm_prev = zeros(1,length(h_lm)); %'topography correction': accounts for chaging shorelines through time.
        delLa_lm_prev = zeros(1,length(h_lm)); %rotational component
        deli_00_prev = 0; %change in ice volume 
        sdelL_lm = zeros(length(ice_time_new)-1,length(h_lm));  %incremental change in Load (ice and ocean mass)
        sdelLa_lm = zeros(length(ice_time_new)-1,length(h_lm)); %incremental change in S (sea surface height) 
        sdelI = zeros(length(ice_time_new)-1,3); %rotational component
        sdelm = zeros(length(ice_time_new)-1,3); %rotational component
        ESL = 0; %eustatic (ice volume equivalent) sea level

        % update new initial topography ('topo_initial' is reset at the end of the loop)
        topo(:,:,1) = topo_initial(:,:,topo_it); 

        % remove the corrected ice model and add predefine (input) ice model back on
        % this needs to be done to recalculate the marine-based ice correction and get a 
        % new 'corrected ice model'
        for i = 1:length(ice_time_new) 
            topo(:,:,i) = topo(:,:,i) - ice_corrected(:,:,i) + ice(:,:,i);
        end

        % do marine-based ice check to recalculate the corrected ice model
        for i = 1:length(ice_time_new)
            if include_ice_check == 'y'
                 % check ice model for floating ice
                 check1 = sign_01(-topo(:,:,i) + ice(:,:,i));
                 check2 = sign_01(+topo(:,:,i) - ice(:,:,i)) .* ...
                     (sign_01(-ice(:,:,i)*rho_ice - (topo(:,:,i) - ice(:,:,i))*rho_water));
                 ice_corrected(:,:,i) = check1.*ice(:,:,i) + check2.*ice(:,:,i);
            else
                % if the floating ice check is set to 'n' then don't change the
                % ice model
                 ice_corrected(:,:,i) = ice(:,:,i);
            end
        end

        % update all topographies with the new corrected ice model
        for i = 1:length(ice_time_new) 
            topo(:,:,i) = topo(:,:,i) - ice(:,:,i) + ice_corrected(:,:,i);
        end

        % assign topography of time 0 and calculate ocean functions for
        % time zero
        topo_0 = topo(:,:,1);
        oc_0 = sign_01(topo_0); %(1=ocean, 0=land). 
        oc0_lm = sphere_har(oc_0,maxdeg,N,P_lm); %convert to spherical harmonics
        ocj_lm_prev = oc0_lm;
        

        % TIME ITERATION
        for t_it = 2:length(ice_time_new) % loop over time

            % Assign topography and ocean function of time t_it to the
            % index j 
            topo_j = topo(:,:,t_it);
            oc_j = sign_01(topo_j);
            ocj_lm = sphere_har(oc_j,maxdeg,N,P_lm); 

            % calculate topography correction (See Kendall et al. 2005 for derivation)
            TO = topo_0.*(oc_j-oc_0);
            TO_lm = sphere_har(TO,maxdeg,N,P_lm); 

            % calculate the change in ice model
            del_ice_corrected = ice_corrected(:,:,t_it) - ice_corrected(:,:,1);
            deli_lm = sphere_har(del_ice_corrected,maxdeg,N,P_lm);
            % calculate the incremental increase in ice volume
            sdeli_00 = deli_lm(1) - deli_00_prev;

            % initial values for convergence
            conv = 'not converged yet';

            % SEA LEVEL EQUATION ITERATION
            for k = 1:k_max  % loop to converge on the change in relative sea level at time t_it

                switch conv

                    case 'converged!' % end loop!

                    case 'not converged yet'

                    % set up initial guess for sea surface height change
                    if k == 1 && topo_it == 1
                        % initial guess of sea level change is just to distribute the
                        % ice over the oceans
                        % (use slightly different initial guess than
                        % Kendall et al. 2005)
                        sdelS_lm(t_it,:) = ocj_lm_prev/ocj_lm_prev(1)*...
                            (-rho_ice/rho_water*sdeli_00 + ...
                            TO_lm(1)-TO_lm_prev(1)) ...
                            - TO_lm - TO_lm_prev;
                    end

                    % add this initial guess to the previous guess from the
                    % last time through the iteration. 
                    delS_lm = delS_lm_prev + sdelS_lm(t_it,:);

                    % calculate total cummulative (from first time step to this one)
                    % change in loading (ice and water)
                    delL_lm = rho_ice*deli_lm + rho_water*delS_lm;
                    % calculate incremental (relative to last time step) change in loading (ice and water)
                    sdelL_lm(t_it-1,:) = delL_lm - delL_lm_prev;

                    
                    % -- calculate viscous deformation contribution -- %

                    % 'beta' contains the viscous love numbers for time t_it,
                    % row index goes over the time increments, column
                    % index goes over lm
                    if t_it == 2 
                        V_lm = zeros(size(T_lm)); %(no viscous contribution at first time step)
                    else
                        for lm_it = 1:length(h_lm) %loop over all love numbers
                            V_lm(lm_it) = beta_l{t_it-1}(:,beta_counter(lm_it))'...
                                * sdelL_lm(1:t_it-2,lm_it); 
                        end
                    end

                    % calculate contribution from rotation
                    if include_rotation == 'y'
                        [delLa_lm, sdelI, sdelm] = calc_rot_visc(delL_lm,...
                            k_el(2),k_el_tide(2),t_it,...
                            beta_konly_l, beta_konly_tide,...
                            sdelI, sdelm);
                        sdelLa_lm(t_it-1,:) = delLa_lm - delLa_lm_prev;
                        if t_it == 2
                            V_lm_T = zeros(size(T_lm));
                        else
                            for lm_it = 1:6 % don't need to loop over all degrees 
                                V_lm_T(lm_it) = beta_tide{t_it-1}(:,beta_counter(lm_it))'...
                                    * sdelLa_lm(1:t_it-2,lm_it);
                            end
                        end   
                        
                        % calculate the sea level perturbation (elastic,
                        % E_lm, and viscous, V_lm)
                        delSLcurl_lm_fl = E_lm .* T_lm .* delL_lm + T_lm .* V_lm + ...
                           1/g*E_lm_T.*delLa_lm + 1/g*V_lm_T;
                    
                    else % if don't include rotation 
                        delSLcurl_lm_fl = E_lm .* T_lm .* delL_lm + ...
                            T_lm .* V_lm;
                    end
                    
                    % convert to spherical harmonics and subtract terms that are part
                    % of the topography to get the 'pure' sea level change
                    delSLcurl_fl = inv_sphere_har(delSLcurl_lm_fl,maxdeg,N,P_lm);
                    delSLcurl = delSLcurl_fl - del_ice_corrected; 


                    % Compute and decompose RO, which represents change in
                    % sea level only over the ocean. (Remember, we ultimately 
                    % calculate 'relative sea level' change everywhere, even on land)
                    RO = delSLcurl.*oc_j;  
                    RO_lm = sphere_har(RO,maxdeg,N,P_lm);

                    % Calculate eustatic sea level perturbation (delta Phi / g)
                    delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*deli_lm(1) ...
                        - RO_lm(1) + TO_lm(1));

                    % Now, we can recalculate a guess for the change in sea surface 
                    % height everywhere in the ocean at this given time. 
                    % (this includes the non-uniform shift in sea surface 
                    % height (RO), the uniform shift (delPhi_g), a correction 
                    % for migrating shorelines (TO) and subtracting all the 
                    % changes in sea surface height from the timesteps before (delS). 
                    sdelS_lm_new = RO_lm + delPhi_g.*ocj_lm - TO_lm ...
                        - delS_lm_prev;

                    % calculate convergence criterion chi
                    chi = abs((sum(abs(sdelS_lm_new)) - sum(abs(sdelS_lm(t_it,:)))) / ...
                        sum(abs(sdelS_lm(t_it,:))) );

                    % check convergence against the value epsilon
                    % If converged, set the variable conv to 'converged!' so that the
                    % calculation exits the loop. If not converged iterate again.
                    if chi < epsilon
                        conv = 'converged!';
                        disp(['Finished time ' num2str(ice_time_new(t_it))...
                        'kyr. Number of iterations ' num2str(k) '. delphi is ' num2str(delPhi_g)])
                       % disp(['Converged after iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    elseif chi < epsilon && k == k_max
                        conv = 'not converged yet';
                        disp(['Finished time ' num2str(ice_time_new(t_it))...
                        'kyr. Run has not converged. Chi is  ' num2str(chi)])
                    else
                        conv = 'not converged yet';
                        %disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
                    end

                    % update incremental change in sea sea surface height
                    % for the next iteration through the sea level loop
                    % (Not necessary if already converged)
                    sdelS_lm(t_it,:) = sdelS_lm_new;

                end

            end % END SLE LOOP

            % * back in time-step loop *
            
            % Set up variables to be used as the 'previous' changes in the
            % next time step. 
            delS_lm_prev = delS_lm;
            TO_lm_prev = TO_lm;
            delL_lm_prev = delL_lm;
            deli_00_prev = deli_lm(1);
            sdelI_00(t_it) = sdeli_00; % record changes through timne for calcaulting GMSL later
            if include_rotation == 'y'
                delLa_lm_prev = delLa_lm;
            end

            % Calculate the ice-equivalent eustatice sea level change here
            % assuming a constant (modern) ocean area
            ESL(t_it) = deli_lm(1)/oc_area * rho_ice/rho_water;

            % Calculate overall perturbation of sea level over oceans
            % (spatially varying field, delSlcurl, and constant offset, delPhi_g)
            delSL = delSLcurl + delPhi_g;

            % Write in topography for next iteration. 
            % (topography is just the negative of sea level)
            topo(:,:,t_it) = - delSL + topo_0;
            ocj_lm_prev = ocj_lm;


        end % END TIME LOOP

        % * back in topography loop *

        % Get the present day topography with ice corrected for
        % marine-based ice.
        topo_pres_ice_corrected = topo_pres - ice(:,:,ind_pres) + ice_corrected(:,:,ind_pres);

        % Compare that modern topography with the present day topography
        % predicted by the GIA simulation.
        topo_diff = max(max(abs(topo(:,:,ind_pres) - topo_pres_ice_corrected)));

        % Determine convergence
        if topo_diff < max_topo_diff
            conv_topo = 'converged!';
            disp(['Converged!! Number of topo iterations ' num2str(topo_it) ...
                '. Topo_diff is ' num2str(topo_diff)])
        else
            conv_topo = 'not converged yet';
            disp(['Not converged. Number of topo iterations ' num2str(topo_it) ...
                '. Topo_diff is ' num2str(topo_diff)])
        end
    end
    
    % Update initial topography for the next time through the topography
    % loop. The new initial topography is set to be the presentday topograhy 
    % minus the change in topography predicted by the GIA simulation. I.e.
    % if the model predicts uplift in Canada, the new initial topography there is
    % corrected downward. (not necessary if already converged). 
    topo_initial(:,:,topo_it+1) = topo_pres_ice_corrected - (topo(:,:,ind_pres) - topo(:,:,1));

end % END OF TOPO LOOP

% calculate relative sea level. This is just the negative of topography, but without
% the change in ice height, so that needs to be subtracted. 
RSL = zeros(size(topo));
for i = 1:length(ice_time_new)
    RSL(:,:,i) = (topo(:,:,ind_pres) - ice_corrected(:,:,ind_pres)) - ...
        (topo(:,:,i) - ice_corrected(:,:,i));
end
toc

% Calculate ESL relative to present. 
ESL = -(ESL - ESL(ind_pres));

%% -------------------------------
% Calculate GMSL 
% --------------------------------
% 
% ESL calculates the ice equivalent sea level change for a fixed ocean
% area. Global mean sea level change (GMSL) calculated here accounts for an 
% ocean area that changes with time following equation 1 of Lambeck et al. 
% (2014). Note that ESL / GMSL are defined differently across the community
% and care should therefore always be taken when using these quantities. 

oc_area_time = zeros(1,length(ice_time_new)); % initialize

%Calculate ocean area through time
for i = 1:length(ice_time_new)
    oc_pres = sign_01(squeeze(topo(:,:,i)));
    oc_area_time(i) = sphere_har(oc_pres,0,N,P_lm);
end
oc_area_rev = flip(oc_area_time);

% Get ice volume changes
sdelI_00_rev = flip(sdelI_00); % sdelI_00 = changes in ice volume over time
sdelI_00_rev = [0 sdelI_00_rev(1:end-1)]; % set it to 0 at time 0

% Calculate GMSL change
GMSL(1) = 0; % at present it's 0
for i = 2:length(ice_time_new)
    temp = rho_ice/rho_water * sdelI_00_rev(i)/oc_area_rev(i); %convert ice volume change to GMSL change
    GMSL(i) = GMSL(i-1) + temp;
end
GMSL = flip(GMSL);


%% -------------------------------
% PLOT RESULTS 
% --------------------------------
%% Plot results: fingerprint

% For fingerprint, look at change in RSL
del_SL = RSL(:,:,2) - RSL(:,:,1);

% calculate the scaling to normalize RSL (it's normalized to be
% one on average, when averaged over the final ocean basin).
FP_normalization_factor = oc0_lm(1)/sphere_har(del_SL.*oc_0,0,N,P_lm); 

% What to plot: (e.g. RSL, ice, topo)
plot_timeslice = squeeze(del_SL*FP_normalization_factor); 

figure
hold on
pcolor(lon_out,lat_out,plot_timeslice)

% plot today's coastlines
lon_shift = [lon_topo(end/2+1:end-1)-360 lon_topo(2:end/2)];
topo_shift = [topo_ice(:,end/2+1:end-1) topo_ice(:,2:end/2)];
contour(lon_shift,lat_topo,topo_shift,[0 0],'w')

shading flat
colorbar 

%% Plot results: time slice RSL

% Set time (in ka)
fig_time = 26; 
ind = find(ice_time_new==fig_time);

% What to plot: (e.g. RSL, ice, topo)
plot_timeslice = squeeze(RSL(:,:,ind)); 

figure
hold on
pcolor(lon_out,lat_out,plot_timeslice)

% plot coastlines at the given time
contour(lon_out,lat_out,topo(:,:,ind),[0 0],'w')

% plot topography (without ice) at the given time
% contour(lon_out,lat_out,topo(:,:,ind)-ice_corrected(:,:,end),[0 0],'w')

% plot present day topography
% contour(lon_out,lat_out,topo(:,:,end),[0 0],'w')

shading flat
colorbar 

%% Plot results: time series 

% Set latitude (-90 to 90) an longitude (-180 to 180)
lon_pt = -53.5;
lat_pt = 66.5;

clear RSL_pt
for i = 1:length(ice_time_new)
    RSL_pt(i) = interp2(lon_out,lat_out,squeeze(RSL(:,:,i)),lon_pt,lat_pt);
    %ice_thk_pt(i) = interp2(lon_out,lat_out,squeeze(ice(:,:,i)),lon_pt,lat_pt);
end

% plot
figure
hold on
plot(ice_time_new,RSL_pt','b')
plot(ice_time_new,ESL','r')
plot(ice_time_new,GMSL','g')
% plot(ice_time_new,ice_thk_pt') % plot thickness of ice sheet at given location
legend('Relative Sea Level','Eustatic Sea Level (ice volume in sea level equivalent)', 'Global Mean Sea Level (as in Lambeck et al. 2014)')
xlim([0 122])
ylabel('Sea Level')
xlabel('Time (ka)')


