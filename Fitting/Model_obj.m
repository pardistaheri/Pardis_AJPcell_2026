function Error = Model_obj(x, p, data)
    p.pest(1:p.NPar) = x; % assign parameters (are in log space)
    
    %% read data
    JO2_Pkd = [data.JO2_Pk(:,1) data.JO2_Pk(:,2) data.JO2_Pk(:,3) data.JO2_Pk(:,4)];
    JO2_Pk_Err = [data.JO2_Pk_Err(:,1) data.JO2_Pk_Err(:,2) data.JO2_Pk_Err(:,3) data.JO2_Pk_Err(:,4)];
    JO2_nmol = [data.JO2_nmol_TC(:,1) data.JO2_nmol_TC(:,2) data.JO2_nmol_TC(:,3) data.JO2_nmol_TC(:,4)];
    dPsi = [data.dPsi(:,1) data.dPsi(:,2) data.dPsi(:,3)];
    RCI_data = [data.RCI(:,1) data.RCI(:,2) data.RCI(:,3) data.RCI(:,4)];

    %% substrates and ADP additions 
    PYR_index = [1 0 0 0];
    GLU_index = [0 1 0 0];
    MAL_index = [1 1 0 0];
    SUC_index = [0 0 1 1];
    ADP_add = [12.5 25 50 100] * 1e-6;
    p.ADPL = length(ADP_add);
    
    options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
    
    % Initialize error
    Error = 0;
    
    % Preallocate cell arrays
    Tc = cell(p.ADPL + 3, p.NSub); % p.Es = 3
    Xc = cell(p.ADPL + 3, p.NSub);
    Jc = cell(p.ADPL + 3, p.NSub);
    Pkc = cell(p.ADPL + 3, p.NSub);
    HPkc = cell(p.ADPL + 3, p.NSub);
    RCIc = cell(p.ADPL + 3, p.NSub);
    
    for i = p.ISub:1:p.NSub 
        X0 = ICs(p); 
        p_tem = p; 
        
        T0 = 0; 
        jj = 1; 
        Es = 3; % extra states 
        cc = 1;
        
        % Initialize Pk array for this substrate
        Pk = zeros(p.ADPL + Es, 1);
        HPk = zeros(p.ADPL + Es, 1);
        RCI = zeros(p.ADPL + Es, 1);
        
        for ii = 1:(p.ADPL + Es)
            %%% Substrate addition 
            if i == 4
                p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
            end 
            
            if ii == (Es - 1)
                X0(p_tem.iMALe) = cc * MAL_index(i) * 5e-3;
                X0(p_tem.iPYRe) = cc * PYR_index(i) * 10e-3;
                X0(p_tem.iGLUe) = cc * GLU_index(i) * 10e-3;
                X0(p_tem.iSUCe) = cc * SUC_index(i) * 7e-3;
            end
            
            %%% ADP addition
            if ii >= Es && ii ~= (p.ADPL + Es)
                X0(p_tem.iADPe) = X0(p.iADPe) + ADP_add(jj);
                jj = jj + 1;
            end
            
            %%% Solving ODEs
            tspan = T0:p.tstep:(T0 + p.time(ii,1));
            [T, X] = ode15s(@ODEs, tspan, X0, options, p_tem);  
            
            T0 = T(end);        
            X0 = X(end,:);
            
            Tc{ii,i} = T;
            Xc{ii,i} = X;
            
            %%% Calculating fluxes - FIXED HERE
            n_fluxes = length(Fluxes(X(1,:), p_tem)); % Get number of fluxes
            J = zeros(length(T), n_fluxes); % Preallocate J with correct dimensions
            
            for zz = 1:length(T)
                flux_result = Fluxes(X(zz,:), p_tem);
                J(zz,:) = flux_result; % This should work now
            end
            
            Jc{ii,i} = J;
            
            s = 14; % CIV flux- OCR
            st = 20; 
            
            %%% finding peaks and length of state 3
            if ii < Es 
                Pk(ii,i) = mean(J(st:end-20, s));
                HPk(ii,i) = mean(J(st:end-20, s));
                RCI(ii,i) = 1;
            else
                Pk(ii,i) = max(J(:, s));
                HPk(ii,i) = 0.5 * Pk(ii,i);
                if ii > 2 % Make sure we have a valid denominator
                    RCI(ii,i) = Pk(ii,i) / Pk(2,i);
                else
                    RCI(ii,i) = 1;
                end
            end
            
            Pkc{ii,i} = Pk(ii,i);
            HPkc{ii,i} = HPk(ii,i);
            RCIc{ii,i} = RCI(ii,i);
        end
         % storing cell variables in vectors 
Jv(:,:,i)=1e9*  [Jc{1,i};   Jc{2,i};    Jc{3,i};    Jc{4,i};    Jc{5,i};    Jc{6,i};]; % nmol/min/mg mito
OCR(:,i)=movmean(Jv(:,14,i),5); % 5points smooth the flux 
        
        % Calculate errors - ensure dimensions match
        valid_indices = 1:min(5, length(Pk)-1);
        
        Err_Pk = sum(((JO2_Pkd(valid_indices,i) - 1e9 * Pk(valid_indices+1,i)) ./ max(JO2_Pkd(:))).^2);
        
        Err_RCI = sum(((RCI_data(valid_indices,i) - RCI(valid_indices+1,i)) ./ max(RCI_data(:))).^2);
        
        min_len = min(length(OCR(:,i)), size(JO2_nmol, 1));
        Err_JO2 = sum(((JO2_nmol(1:min_len,i) - OCR(1:min_len,i)) ./ max(JO2_nmol(:))).^2);
        
        Error = Error + Err_JO2 + Err_Pk + Err_RCI;
    end
end