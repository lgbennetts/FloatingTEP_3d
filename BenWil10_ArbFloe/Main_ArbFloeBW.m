% function Main_ArbFloeBW
%
% INPUT:
%
% fxn_handle  = @CrookedEgg or @EllipticalCurve;
% fraxn       = parameter for shape function (usually = 1)
% VM_vec      = vector of Vert_Dim
%               Vert_Modes = The number of natural vertical modes used (>0)
% GM_vec0     = [GM_vec, Gegenbauer yes/no]
%               GM_vec = vector of DimG
%               DimG = The number of Gegenbauer Modes used (>0)
%               Gegenbauer yes/no = use the Gegenbauer polys at vertical interface (0=no)
% SM_vec      = vector of s_Modes
%               s_Modes = The limt of the Fourier expansion in horizontal plane (any nat num)
%               total number = 2*s_Modes+1
% forceinput  = {forcing type, value}
%               forcing type = 'p' wave period; 'k' wavenumber water; 'i' wavenumber ice
%               inc_angle = angle incident wave makes with x-axis
% radius      = radius of floe (>0)
% ellipse     = 2x1 vector of stretch in x/y directions
% rotation    = rotate floe (in degs)
% translation = translate floe from origin 2x1 vector
% fig         = specify figure (if plotting)
%
% OUTPUT:
%
% E0     = scattered energy
% k_out  = [water wavenumber, ice wavenumber]
% FF_Amp = far-fied amplitude (fn of theta)
% th_vec = abscissa for FF_Amp
% mx     = maximum average bending (laplace*w)
%
% MISC:
%
% wt_typ   = weighting on vertical modes
%            'reg'  -> 1/cosh(k*h)
%            'norm' -> 1/<w_i>
%            'none' -> 1
% Nint     = mesh points
% Norm_fac = normalisation of lengths
%
% L Bennetts 2010 / Dunedin
%
% Revision history:
%                   L Bennetts Nov 2013

function  ...  % [E0, k_out, FF_Amp, th_vec, mx] =
 Main_ArbFloeBW(fxn_handle, fraxn, VM_vec, GM_vec0, SM_vec,...
 forceinput, inc_angle, radius, ellipse, rotation, translation, fig)

if ~exist('fxn_handle','var'); fxn_handle=@SuperCircle; end
%if ~exist('fxn_handle','var'); fxn_handle=@EllipticalCurve; end
if ~exist('fraxn','var'); fraxn=4; end
if ~exist('VM_vec','var'); VM_vec=1; end
if ~exist('GM_vec0','var'); GM_vec0=[1,0]; end
if ~exist('SM_vec','var'); SM_vec=2^4; end
if ~exist('forceinput','var'); forceinput={'p',5}; end
if ~exist('inc_angle','var'); inc_angle=0; end
if ~exist('radius','var'); radius=5; end
if ~exist('ellipse','var'); ellipse=[1,1]; end
if ~exist('rotation','var'); rotation=0; end
if ~exist('translation','var'); translation=[0,0]; end
if ~exist('fig','var'); fig=0; end

if ~exist('DO_COMM','var'); DO_COMM=1; end

wt_typ   = 'reg';
Nint     = 4*300;
Norm_fac = 1;

%%

Tol_vec(1) = 1e-16; % Real root error %
Tol_vec(2) = 1e-16; % Imag root error (init) %
Tol_vec(3) = 1e-4;  % Tol on const soln %

%% Define geom parameters

parameter_vector = Param_vec(0); % (just to get started)

thickness = 1; bed = 20;
draught = parameter_vector(2)*thickness/parameter_vector(1);

%% Define inc wavelength

if forceinput{1}=='p'
 freq = 2*pi/forceinput{2};
elseif forceinput{1}=='k'
 freq = FindFreq_FS(parameter_vector, forceinput{2}, bed);
elseif forceinput{1}=='i'
 freq = FindFreq_ice(forceinput{2}, ...
  (parameter_vector(2)/parameter_vector(1))*thickness, ...
  (parameter_vector(4))*(thickness^3)/...
  (12*parameter_vector(1)*(1-(parameter_vector(3)^2))*parameter_vector(5)), bed-draught);
end

%% Normalise

thickness = thickness/Norm_fac; draught = draught/Norm_fac;
%lam = lam/Norm_fac;
freq = freq/sqrt(Norm_fac);
bed = bed/Norm_fac;
radius = radius/Norm_fac;

al = (parameter_vector(2)/parameter_vector(1))*thickness; % - mass
be = (parameter_vector(4)/Norm_fac)*(thickness^3)/...
 (12*parameter_vector(1)*(1-(parameter_vector(3)^2))*parameter_vector(5)); % - flex rigid

% freq = FindFreq_FS(parameter_vector, 2*pi/lam, bed); % - non-dim frequency
parameter_vector = Param_vec(freq); % - densities, etc
parameter_vector(4) = parameter_vector(4)/Norm_fac;

%display(['freq=', num2str(freq), ', period=', num2str(2*pi/freq)])

parameter_vector = [parameter_vector(1:6), bed, draught, thickness, bed-draught, al, be];

GM_vec = GM_vec0(1:end-1); Geg_yn = GM_vec0(end); clear GM_vec0

% display([al, be])

%% Define truncations

% figure(fig); hold on
count=0;

for loop_vm=1:length(VM_vec)
 Vert_Modes = VM_vec(loop_vm);
 for loop_g=1:length(GM_vec)
  DimG = GM_vec(loop_g);
  for loop_s=1:length(SM_vec)
   s_Modes = SM_vec(loop_s);
   count = count+1; %col = col_vec{count};
   
   % Vert_Modes = 1; % 1 is minimum
   % s_Modes    = 1; % 0 is minimum, and adds in Fourier +/- pairs
   
   %% Define matricies
   
   % wvnos in ice (kk) and water (k0)
   kk = zeros(Vert_Modes,1); k0 = zeros(Vert_Modes,1);
   % weight (normalising fn) attached to each vert mode
   wt = zeros(Vert_Modes,1); wt_0 = zeros(Vert_Modes,1);
   
   for loop_Dim = 1:Vert_Modes
    
    kk(loop_Dim) = GetRootsMMA_PWC(parameter_vector, loop_Dim, Tol_vec);
    wt(loop_Dim) = weight_PWC(parameter_vector, kk(loop_Dim), wt_typ);
    
    k0(loop_Dim) = GetRootsMMA_FS_PWC(parameter_vector, loop_Dim, Tol_vec);
    wt_0(loop_Dim) = weight_0_PWC(parameter_vector, k0(loop_Dim), wt_typ);
    
   end
   
   % apx complex roots (Dimension dependent)
   [mu_0, mu_1] = mu_new_PWC(parameter_vector, Vert_Modes, kk, wt);
   
   mat_A = matrix_A_PWC(parameter_vector, Vert_Modes, kk, wt, parameter_vector(10)); % in ice
   mat_A0 = matrix_A_PWC(parameter_vector, Vert_Modes, k0, wt_0, parameter_vector(7)); % in water
   mat_A0 = diag(mat_A0);
   
   C_mat = mat_C_dr_PWC(parameter_vector, length(kk), kk, mu_0, mu_1,...
    wt, parameter_vector(12));
   % matrix that controls motion in ice
   % the `complicated' weight stuff is all done in C_mat
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Vertical modes %%%
   if ~Geg_yn
    [mat_W0, mat_W] = jump_W_PWC(parameter_vector, Vert_Modes, DimG, kk,...
     k0, kk, wt_0, wt);
    %%% Gegenbauer %%%
   else
    if draught==0
     [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, bed, draught, 0.5, DimG);
    else
     [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, bed, draught, 1/6, DimG);
    end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %% - Incident wave - %
   
   % inc_angle = 0; % - (wrt x-axis) keep this FIXED as can just rotate the pool
   inc_amp = parameter_vector(6)/wt(1)/kk(1)/sinh(kk(1)*(bed-draught))...
    /Norm_fac; % 1; % - FIXED to normalise to 1m amplitude (in ice)
   
   %% - Solve integral eqn - %
   
   matrices={mat_W0,mat_W,C_mat,mat_A};
   vectors={k0,[kk; mu_0; mu_1],mat_A0};
   parameters={al,be,parameter_vector(6),parameter_vector(3), inc_angle};
   % fxn_handle = @CrookedEgg;
   % fxn_handle = @EllipticalCurve;
   % fraxn=1;
   scaler=radius*ellipse;
   % rotation=0;
   % translation=[0 0];
   srt={scaler,rotation,translation};
   edge_vars={fxn_handle,{fraxn},srt};
   %Nterms = 5;
   ff_Modes = s_Modes;
   
   if DO_COMM
    cprintf(0.4*[1,1,1],['>>> ' char(fxn_handle) '; radius=' num2str(radius) ...
     '; stretch=' num2str(ellipse(1)) ',' num2str(ellipse(2)) '\n'])
    cprintf(0.4*[1,1,1],['>>> wavelength=' num2str(2*pi/k0(1)) ': ' ...
     num2str(inc_angle) ' degs\n'])
    cprintf(0.4*[1,1,1],['>>> thickness=' num2str(thickness)])
    cprintf(0.4*[1,1,1],['>>> Modes: vert=' int2str(Vert_Modes) ...
     '; ang=' int2str(s_Modes) '\n'])
    cprintf(0.4*[1,1,1],['>>> Integration points ' int2str(Nint) '\n'])
   end
   
   [dum_Rp,dum,LC,Phi,Phi_dn,Psi,Psi_dn,svec]=...
    FLOE_sys_solve_smooth_vG(matrices,vectors,...
    parameters,edge_vars,Vert_Modes-1,DimG-1,s_Modes,ff_Modes,Nint); %,...
   %fig,col);
   
   clear dum
   
   % for loop_chi=-s_Modes:s_Modes
   %     chi_mat(s_Modes+1+loop_chi,:) = exp(1i*loop_chi*pi*svec/LC);
   % end
   % chi_mat = chi_mat/2/LC;
   %
   % wnew = real(Psi(Vert_Modes+1,:)*chi_mat);
   %
   % if count~=1
   % disp(L2reler(wold,wnew))
   % end
   %
   % wold = wnew;
   %
   % if 0
   %  plot(svec/LC, wnew, col);
   % end
   
   %% - Checks - %
   
   %% -- Far-field energy -- %%
   
   % Energy = inc_amp*mat_A(1,1)*(freq/sqrt(Norm_fac))*sum(abs(Rp).^2)
   
   % Rp = dum_Rp(1,:).*exp([-ff_Modes:ff_Modes]*pi/2i);
   % E0 = sum(abs(Rp).^2); E1 = real(sum(Rp));
   
   inc_ang = inc_angle*pi/180;
   var_th = GrafThetaVal(cos(inc_ang), sin(inc_ang), [1, 1]);
   Rp = dum_Rp(1,:); %inc_amp*wt(1)*
   Ip = exp(1i*[-ff_Modes:ff_Modes]*var_th); %inc_amp*wt(1)*
   E0 = sum(abs(Rp).^2); E1 = real(sum(conj(Ip).*Rp));
   
   %% - Far-field amplitude
   
   th_vec = linspace(-pi,pi,501);
   Rp = Rp.*exp([-ff_Modes:ff_Modes]*pi/2i);
   FF_Amp = Rp*exp(1i*[-ff_Modes:ff_Modes].'*th_vec);
   
   %% - check - %
   
   if DO_COMM
    if abs(E0-E1)>1e-3
     cprintf('magenta',['energy check ' num2str(abs(E0-E1)) '\n']);
    end
   end
   
   E0 = E0/2/pi;
   
   %% - For the scattering X sect - %
   
   % if count~=1
   %  disp(abs(Eold-E0)/E0)
   % end
   %
   % Eold = E0;
   
   %% - For the amplitude plots - %
   
   Fnew = abs(FF_Amp);
   
   if count~=1
    disp(L2reler(Fold,Fnew))
   end
   
   Fold = Fnew;
   
  end % - s_Modes loop
 end % - GDim loop
end % - Vert_Dim loop

k_out(1) = k0(1); k_out(2) = kk(1);

Psi = inc_amp*Psi;
Psi_dn = inc_amp*Psi_dn;

if exist('fig','var')
 if fig
  mx = Max_Floe_new(edge_vars, 2*25, 2*20, Psi, Psi_dn,...
  [kk.', mu_0, mu_1], C_mat, Vert_Modes, s_Modes, 50, fig);
 
  mx = mx/be;
 else
  mx=0;
 end
else
 mx=0;
end

if DO_COMM
 cprintf('blue',['>> ... scattered energy=' num2str(E0) '\n'])
end

return

%% -- disp in pool
Phi = inc_amp*Phi;
Phi_dn = inc_amp*Phi_dn;
Psi = inc_amp*Psi;
Psi_dn = inc_amp*Psi_dn;

% - Contour plots:

Plot_Floe(edge_vars, 50, [2*10,2*25], Psi, Psi_dn, Phi, Phi_dn, [kk.', mu_0, mu_1],...
 C_mat, k0, Vert_Modes, s_Modes, [inc_amp, inc_ang], [{bed},{wt_0}], fig)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameter_vector = Param_vec(freq)

liquid_density = 1025; %1000; %  %in kg/m^3
ice_density = 922.5; %900; %  %in kg/m^3 %
poissons_ratio = 0.3; %1/3; %
youngs_mod = 6e9; %5000000000;

gravity = 9.81; %1; % in m/sec^2
kappa = (freq^2)/gravity;

parameter_vector = [liquid_density, ice_density, poissons_ratio, ...
 youngs_mod, gravity, kappa];

return

% function that produces frequency from given waveno in ice-covered region

function freq = FindFreq_ice(kk, al, be, HH)

gravity = Param_vec(0);

freq = (1+be*(kk^4))*kk.*tanh(kk.*HH)/(1+al*kk.*tanh(kk.*HH));

freq = sqrt(gravity(5)*freq);

return

% function that produces frequency from given waveno in free-surface region

function freq = FindFreq_FS(parameter_vector, kk, HH)

freq = kk.*tanh(kk.*HH);

freq = sqrt(parameter_vector(5)*freq);

return

%

function w = weight_PWC(parameter_vector, root, typ)

if strcmp(typ,'norm')
 w = normalising_weight_PWC(parameter_vector, root);
elseif strcmp(typ,'none')
 w = 1;
elseif strcmp(typ,'reg')
 H = parameter_vector(10);
 w =  sech(root*H);
end

return

%

function w0 = weight_0_PWC(parameter_vector,  root, typ)

if strcmp(typ,'norm')
 w0 =  normalising_weight_0_PWC(parameter_vector, root);
elseif strcmp(typ,'none')
 w0 = 1;
elseif strcmp(typ,'reg')
 h = parameter_vector(7);
 w0 =  sech(root*h);
end

return