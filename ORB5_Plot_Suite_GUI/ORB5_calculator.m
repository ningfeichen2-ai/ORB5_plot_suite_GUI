function plasma = ORB5_calculator(Z,mu,filename)
% GET_PLASMA_PARAMETERS Calculates CGS plasma units and ORB5 scaling factors.
%
% Inputs:
%   lx    - 2a/rho_i (from ORB5 input)
%   Btor  - Magnetic field B0 [Tesla]
%   beta  - Plasma beta (ORB5 definition: half of usual beta)
%   amid  - Minor radius [m]
%   r0mid - Major radius [m]
%   Z     - Charge number
%   mu    - Mass number
%   q0    - Safety factor
%
% Output:
%   plasma - A struct containing CGS parameters and transformation coefficients

    if nargin < 3, filename = fullfile(pwd,'input'); end
    validateattributes(Z,{'numeric'},{'scalar','positive','finite'});
    validateattributes(mu,{'numeric'},{'scalar','positive','finite'});
    %% 1. Physical Constants (CGS)
    e_cgs = 4.8032e-10;         % Elementary charge [statC]
    c_cgs = 3e10;               % Speed of light [cm/s]
    k_ev_to_erg = 1.602e-12;    % erg/eV conversion factor

    %% 2. Basic Geometry and Field Conversion
    plasma.a_cgs  = get_orb5_val('equil','a_mid',filename) * 10^2;    % m -> cm
    plasma.R0_cgs = get_orb5_val('equil','r0_mid',filename)* 10^2;   % m -> cm
    plasma.B0_cgs = get_orb5_val('equil','btor0',filename) * 10^4;    % Tesla -> Gauss
    plasma.lx = get_orb5_val('equil','lx',filename);
    plasma.beta = get_orb5_val('equil','beta',filename);


    %% 3. Derived Plasma Parameters (CGS)
    % Sound Larmor radius [cm]
    plasma.rhos = 2 * plasma.a_cgs / plasma.lx; 
    
    % Ion cyclotron frequency [rad/s]
    plasma.wci = 9.58e3 * Z * (mu^(-1)) * plasma.B0_cgs; 
    
    % Temperature [eV]
    plasma.Te = 4 * plasma.a_cgs^2 * plasma.B0_cgs^2 * Z^2 / (mu * 10^4 * 1.02^2 * plasma.lx^2); 
    
    % Sound velocity [cm/s]
    plasma.Cs = 9.79e5 * sqrt(Z * plasma.Te / mu); 
    
    % Equilibrium number density [cm^-3]
    plasma.n0 = 4.96e10 * plasma.beta * plasma.B0_cgs^2 / plasma.Te; 
    
    % Alfven velocity [cm/s] and frequency [rad/s]
    plasma.vA = 2.18e11 * (mu^-0.5) * (plasma.n0^-0.5) * plasma.B0_cgs; 
    plasma.wA0 = plasma.vA / plasma.R0_cgs;
    plasma.wCs = plasma.Cs / plasma.R0_cgs;
    plasma.wci_wA0 = plasma.wci/plasma.wA0;
    plasma.wci_Cs = plasma.wci/(plasma.Cs/plasma.R0_cgs);

    %% 4. Transformation Coefficients (ORB5 -> SI/CGS)
    
    % Vector Potential (Apar) transforms
    % Coefficient for e*deltaA/Te
    plasma.coeff_A_to_eTe = 10^6 * (plasma.rhos * plasma.B0_cgs / 10^6) * (e_cgs / (k_ev_to_erg * plasma.Te)); 
    
    % Coefficient for e*psi/Te (for AEs)
    plasma.coeff_A_to_eTepsi = plasma.coeff_A_to_eTe * plasma.vA / c_cgs; 
    
    % SI and CGS units for Vector Potential
    plasma.coeff_A_to_SI  = (plasma.rhos * plasma.B0_cgs / 10^6); 
    plasma.coeff_A_to_CGS = plasma.coeff_A_to_SI * 10^6;

    % Scalar Potential (Phi) transforms
    plasma.coeff_Phi_to_SI  = plasma.Te; 
    plasma.coeff_Phi_to_CGS = (k_ev_to_erg * plasma.Te) / e_cgs;
    
  
end
