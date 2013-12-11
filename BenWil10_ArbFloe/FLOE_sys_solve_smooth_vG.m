function [Rp,Edge,LC,Phi,Phi_n,Psi_hat,Psi_n_hat,svec]=...
             FLOE_sys_solve_smooth_vG(matrices,vectors,...
               parameters,edge_vars,N,NG,Nterms,Nmodes_out,Nint, fig, col)
	       	       
%% CALL: FLOE_sys_solve_smooth(matrices,vectors,parameters,...
%%        edge_vars,Nterms)
%% matrices={Pminus,Pplus,CMAT,Aplus}
%% vectors={k_in,k_out,dAminus}
%% parameters={alf,beta,sigma,nu,theta_inc in degrees}
%% edge_vars={fxn_handle,fxn_parameters,srt}
%% srt is best seen by looking below in the test inputs

% - 20.08.09 - rewritten by LB for floe prb - %
% - adapted from POOL_sys_solve_smooth_vG - %
% - 15.07.09 - rewritten by LB for Gegenbauer poly exp at interface - %

USE_EXP=1;%% use exponentials for basis,
          %% else use cos's & sin's.

if nargin==0%% CREATE SOME TEST INPUTS:
  N=0;
  Nterms=4;
  %%
  Pminus=eye(N+1);
  Pplus=eye(N+1);
  CMAT=eye(N+3);
  Aplus=eye(N+1);
  matrices={Pminus,Pplus,CMAT,Aplus};
  %%
  dAminus=diag(Aplus);
  k_in=1.1*[1; (1:N)'*i];
  k_out=2*[k_in;1+i;-1+i];
  vectors={k_in,k_out,dAminus};
  %%
  alf=1;
  beta=1;
  sigma=.5;
  nu=.3;
  th_inc=0;
  parameters={alf,beta,sigma,nu,th_inc};
  %%
  edge_fxn=@CRKprof_circarc;%% creates an elliptical arc shape
  fraxn=1;%% full ellipse
  radius=10;
  scaler=radius*[1 1];%% in general this is [x-axis,y-axis]
  rotation=0;%% can also rotate it
  translation=[0 0];%% can also translate it
  srt={scaler,rotation,translation};
  edge_vars={edge_fxn,{fraxn},srt};
  %%
  Nmodes_out=5;
end

alf=parameters{1};
bet=parameters{2};
sigma=parameters{3};
nu=parameters{4};
theta_inc=parameters{5}*pi/180;
bet0=bet*(1-nu);

%% GET QUADRATURE POINTS, WEIGHTS ETC:
Nunc=2*Nterms+1;%% using cosines and sines in Galerkin scheme.
u_coeffs=zeros(Nunc,N+1);
v_coeffs=u_coeffs;
w_coeffs=u_coeffs(:,1);
dnw_coeffs=w_coeffs;
xtra_outs=[];
%%
%Nint=300;
[tt,ww]=OP_numint_legendre(Nint);%% use Gauss-legendre points
                          %% and weights as everything smooth.
[xyvecs,ds_dt,th_vec,dtheta_ds,d2s_dt2,d2theta_ds2,...
   d2xy_ds2,LC,svec]=CRK_get_rsdtheta(edge_vars,tt);

if 0 % - WHEN LOOKING AT CORNERS!!!!!!!!!!!!
disp('################################');
disp(['corner angle=',num2str(abs(th_vec(1)-th_vec(Nint)-pi)/pi),'*pi']);
disp('################################');
end

xvec=xyvecs(1,:).';
yvec=xyvecs(2,:).';
[X0,X]=meshgrid(xvec,xvec);
[Y0,Y]=meshgrid(yvec,yvec);
[R,varTH]=GEN_polar_coords(X-X0,Y-Y0);
[TH0,TH]=meshgrid(th_vec,th_vec);
jnz=find(R);

if 1%% CHECK SHAPE & TANGENT/NORMAL VECTORS:
  %dtheta_ds,1/radius
  figure(fn_getfig)
  plot(xvec,yvec), hold on;
  svectors=[cos(th_vec),sin(th_vec)].';
  nvectors=[sin(th_vec),-cos(th_vec)].';
  for j=1:10:Nint
    svector_x=[0;svectors(1,j)];
    svector_y=[0;svectors(2,j)];
    plot(xvec(j)+svector_x,yvec(j)+svector_y,'m');
    %%
    nvector_x=[0;nvectors(1,j)];
    nvector_y=[0;nvectors(2,j)];
    plot(xvec(j)+nvector_x,yvec(j)+nvector_y,'g');
  end
  hold off;
  daspect([1 1 1])
  drawnow
%   return;
end

%% GET INNER PRODUCT MATRICES:
%ds_dt,pause
ww=ww.*ds_dt/LC;
if USE_EXP%% exponentials:
  [IP,hn,Wn]=GEN_inprod_exp(svec/LC,ww,Nterms); %_vOld
  hn=hn*LC;
  %%
  nn_unc=(-Nterms:Nterms).';
  kunc=pi/LC*nn_unc;
  Ddiff=diag(i*kunc);
else%% cos's & sin's:
  [IP,hn,Wn]=GEN_inprod_cos_sin(svec/LC,ww,Nterms);
  hn=hn*LC;
  %%
  nn_unc=1:Nterms;
  kunc0=pi/LC*nn_unc';
  Ddiff(nn_unc+1,nn_unc+1+Nterms)=-diag(kunc0);
  Ddiff(nn_unc+1+Nterms,nn_unc+1)=diag(kunc0);
  %%
  nn_unc=[0;nn_unc';nn_unc'];
  kunc=pi/LC*nn_unc;
end
%hn=[2*LC;LC+0*(1:Nunc-1)'];

%% CALC HANKEL FUNCTIONS & THEIR NORMAL DERIVATIVES,
%% THEN TAKE THE INNER PRODUCTS WITH THEM:
k_in=vectors{1};
k_out=vectors{2};
N=length(k_in)-1;
matrices4U={R,varTH,TH};
vectors4U={k_in,k_out,d2xy_ds2,svec/LC};

%% need to add 1/2\pi*log(R) singularity into U's
%% ( this is the same for all the k's ):
[S0,SS]=meshgrid(svec,svec);
if USE_EXP
  Exp_term=abs( 1-exp(pi*i/LC*(SS-S0)) );
  %% coefficients for 1/(2*pi)*log[1-exp(i*(s-s_0)/LC)]:
  ln_coeff=-1/4/pi./(1:Nterms).';
  ln_coeff=[flipud(ln_coeff);0;ln_coeff];
else
  Exp_term=abs( 1-exp(pi*i/LC*(SS-S0)) );
  %% coefficients for 1/(2*pi)*log|1-exp(i*(s-s_0)/LC)|:
  ln_coeff=-1/2/pi./(1:Nterms).';
  ln_coeff=[0;ln_coeff;ln_coeff];
end
log_corr=diag( 0*ds_dt - 1/2/pi*log(pi/LC) );
log_corr(jnz)=1/2/pi*log( R(jnz)./Exp_term(jnz) );
%1/2/pi*log(radius),return
%%
log_corr_mn=diag(ln_coeff)+IP*log_corr*IP';

if 0%% CHECK THE BOUNDED PART OF
    %% THE LOG CORRECTION TERM'S EXPANSION:
  nchk=GEN_rand(Nint)
  nnfm=(-Nterms:Nterms)';
  Am=IP*log_corr*IP'*exp(-i*nnfm*pi*svec(nchk)/LC);
  plot(svec/LC,log_corr(:,nchk)), hold on;
  Lap=GEN_interp_exp(svec/LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot(svec/LC,Lap,'--r'), hold off;
  return;
elseif 0%% CHECK THE UNBOUNDED PART OF
        %% THE LOG CORRECTION TERM'S EXPANSION:
  nchk=GEN_rand(Nint)
  nnfm=(-Nterms:Nterms)';
  Am=diag(ln_coeff)*exp(-i*nnfm*pi*svec(nchk)/LC);
  plot(svec/LC,1/2/pi*log(Exp_term(:,nchk))), hold on;
  Lap=GEN_interp_exp(svec/LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot(svec/LC,Lap,'--r'), hold off;
  return;
end

%%
ipstuff={IP,log_corr_mn};
[U_in_mn,U_out_mn,dn0U_in_mn,dn0U_out_mn]=...
POOL_hankel_derivs2(matrices4U,vectors4U,ipstuff);

%% Solve the deriv JCs %%

mat_W0 = matrices{1};
mat_W = matrices{2};
mat_A0 = vectors{3};
mat_A = matrices{4};

mat_JC_ds0 = diag(1./mat_A0)*mat_W0;
mat_JC_ds = mat_A\mat_W;

%% CALC KERNEL MATRIX FROM INNER PROBLEM:
% - assuming that the number of modes of Fourier modes used in Phi and u
% are the same

Jvec=(1:Nunc)';
jv_w=Jvec+(N+1)*Nunc;
jv_del2w=jv_w+Nunc;
jh_w=Jvec+2*(N+1)*Nunc;
jh_dnw=jh_w+Nunc;

IE_Mat_Phi = zeros((N+1)*Nunc);
IE_Mat_u0 = zeros((N+1)*Nunc,(NG+1)*Nunc);

% - FORCING TERMS:
forcing_in=IE_Mat_Phi(:,1);%% this is zero on the inside,
                        %% so leave as is.
kap0=k_in(1);
ff=exp(1i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
fm=IP*ff;
% ------------------- %
                       
Mdiag=.5./hn;%% NB factor 1/2;
   %% this matrix goes on the diagonal of the main matrix
   %% - hn is because we take norm of chi_n so fac of 2LC (LB 11.11.08)
   
% - construct IE_Mat_Phi*Phi = forcing_in + IE_Mat_u0*U - %

v1 = 1:N+1; 
for loop_s1=1:Nunc
forcing_in(v1) = [fm(loop_s1); zeros(N,1)];
v2 = 1:N+1; u2 = 1:NG+1;
for loop_s2=1:Nunc
    IE_Mat_Phi(v1,v2) = diag(squeeze(dn0U_in_mn(loop_s1,loop_s2,:)));
    IE_Mat_u0(v1,u2) = diag(squeeze(U_in_mn(loop_s1,loop_s2,:)))*mat_JC_ds0;
    v2=v2+N+1; u2=u2+NG+1;
end
IE_Mat_Phi(v1,v1) = IE_Mat_Phi(v1,v1) + Mdiag(loop_s1)*eye(N+1);
v1=v1+N+1; 
end

%% CALC KERNEL MATRIX FROM OUTER PROBLEM:
CMAT=matrices{3};
invCMAT=inv(CMAT);
%%
IE_Mat_tilPsi = zeros((N+3)*Nunc);
IE_Mat_u = zeros((N+3)*Nunc,(NG+1)*Nunc);

% - construct IE_Mat_Psi*Psi = forcing_out + IE_Mat_u*U - %

Id = eye(N+3);
tilI = Id(:,1:N+1);
tilIw = Id(:,N+2:N+3); clear Id

% - FORCING TERMS:
forcing_out=IE_Mat_tilPsi(:,1);
% kap0=k_out(1);
% ff=exp(i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
% fm=IP*ff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LB (03.12.08) for debugging only!!!! %%%%%%%%
% ff_dx=i*kap0*cos(theta_inc)*exp(i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
% ff_dy=i*kap0*sin(theta_inc)*exp(i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
% ff_dn=diag( sin(th_vec) )*ff_dx - diag( cos(th_vec) )*ff_dy;
% fm_dn=conj(IP)*ff_dn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0%% TEST FORCING TERM FOR A CIRCLE:
  nnfm=(-Nterms:Nterms)';
  theta_inc%% NB need to change expansion if theta_inc~=0
  fm2=i.^nnfm.*besselj(nnfm,kap0*radius);
  tst=[fm,fm2]
  tst(:,1)./tst(:,2)
  %%
  plot(tt,[real(ff) imag(ff)]), hold on;
  fap=GEN_interp_exp(svec/LC,fm);
  plot(tt,[real(fap) imag(fap)],'--r');
  fap2=GEN_interp_exp(svec/LC,fm2);
  plot(tt,[real(fap2) imag(fap2)],'--g'), hold off;
  return;
end
%  forcing_out(Jvec)=fm;
%  amp_disp=1;%% NB check this
%  forcing_out(Jvec+(N+1)*Nunc)=sigma*amp_disp*fm;
%  forcing_out(Jvec+(N+2)*Nunc)=-k_out(1)^2*sigma*amp_disp*fm;

%% BACK TO KERNEL MATRIX:

%% matrices to apply bending moment condition:
%% del2w_coeffs=MB1*w_coeffs +MB2*dnw_coeffs 
MB1=diag(-bet0*kunc.^2);%% from \beta_0*\pa_s^2w
if 1
  %MB2=diag(hn)*IP*diag(bet0*dtheta_ds)*Wn*diag(1./hn);
  %%%%%%% LB (27.01.08) %%%%%%%%
  MB2=diag(hn)*IP*diag(bet0*dtheta_ds)*Wn*diag(1./hn);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% from bet0*dtheta_ds*\pa_nw
else
  if USE_EXP
    bhat=GEN_inprod_exp_vOld(svec/LC,ww,2*Nterms,dtheta_ds);
    JJ=0:2*Nterms;
    for j=0:2*Nterms
      B(j+1,:)=bhat(JJ+2*Nterms-j+1).';
    end
  else
    disp('make B matrix for cosines and sines')
  end
  MB2=bet0*B;
  %% from bet0*dtheta_ds*\pa_nw
end

%% matrices to apply shearing condition:
%% dn_del2w_coeffs=MS1*w_coeffs +MS2*dnw_coeffs
MS1=Ddiff*MB2*Ddiff;%% from \beta_0\pa_s*dtheta_ds*\pa_sw
%  %%%%%%% LB (12.11.08) %%%%%%%%
%  MS1=-Ddiff*MB2*Ddiff;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
MS2=-MB1;%% from -\beta_0*\pa_s^2\pa_nw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Psi_inds = zeros(Nunc*(N+1),1);
w_inds=zeros(Nunc,1); wn_inds=w_inds;

v1 = 1:N+3; 
for loop_s1=1:Nunc

%     forcing_out(v1) = CMAT*[fm(loop_s1); zeros(N+2,1)];

    v2 = 1:N+1; u2 = 1:NG+1; w2 = N+2:N+3;

    Psi_inds([1:N+1]+(loop_s1-1)*(N+1)) = v1(1):v1(N+1); % - use this later
    w_inds(loop_s1) = v1(N+2); wn_inds(loop_s1) = v1(N+3);

    for loop_s2=1:Nunc
        if loop_s1==loop_s2
            dnm=1; else dnm=0;
        end
        mat_BS = [dnm, 0; MB1(loop_s1,loop_s2), MB2(loop_s1,loop_s2)];
        mat_BS_dn = [0, dnm; MS1(loop_s1,loop_s2), MS2(loop_s1,loop_s2)];

        M0 = CMAT*diag(squeeze(dn0U_out_mn(loop_s1,loop_s2,:)))*invCMAT;
        M1 = CMAT*diag(squeeze(U_out_mn(loop_s1,loop_s2,:)))*invCMAT;

        IE_Mat_tilPsi(v1,v2) = -M0*tilI;
        IE_Mat_u(v1,u2) = -M1*tilI*mat_JC_ds;

        IE_Mat_tilPsi(v1,w2) = IE_Mat_tilPsi(v1,w2) + Mdiag(loop_s1)*tilIw*mat_BS;
        
        w3 = N+2:N+3;

        for loop_s3=1:Nunc
            
            if loop_s2==loop_s3
                dnm=1; else dnm=0;
            end
            mat_BS = [dnm, 0; MB1(loop_s2,loop_s3), MB2(loop_s2,loop_s3)];
            mat_BS_dn = [0, dnm; MS1(loop_s2,loop_s3), MS2(loop_s2,loop_s3)];

            IE_Mat_tilPsi(v1,w3) = IE_Mat_tilPsi(v1,w3) + ...
                -M0*tilIw*mat_BS + M1*tilIw*mat_BS_dn;
            
            w3=w3+N+3;

        end

        v2=v2+N+3; u2=u2+NG+1; w2=w2+N+3;
    end
    IE_Mat_tilPsi(v1,v1(1):v1(N+1)) = IE_Mat_tilPsi(v1,v1(1):v1(N+1)) + Mdiag(loop_s1)*tilI;
    v1=v1+N+3;
end

%% - Solve for Phi, Psi in terms of u (and inc wave) - %%

% - size Nunc*(N+1) x Nunc*(NG+1)

Phi_u0 = IE_Mat_Phi\IE_Mat_u0;
tilPsi_u = IE_Mat_tilPsi\IE_Mat_u;
Psi_u = tilPsi_u(Psi_inds,:);

% - size Nunc*(N+1) x 1

Phi_I0 = IE_Mat_Phi\forcing_in;

% - size Nunc*(N+1) x 1

tilPsi_I = IE_Mat_tilPsi\forcing_out;
Psi_I = tilPsi_I(Psi_inds,:);

% - also define...

% - size 1 x Nunc*(NG+1)
w_u = tilPsi_u(w_inds,:);
wn_u = tilPsi_u(wn_inds,:);
% - size 1 x Nunc*(N+3)
w_I = tilPsi_I(w_inds,:);
wn_I = tilPsi_I(wn_inds,:);

% - Solve for u - %

Big_W0T = zeros(Nunc*(NG+1),Nunc*(N+1)); Big_WT=Big_W0T;
Big_JCds0 = Big_W0T'; Big_JCds=Big_JCds0;

v1 = 1:N+1; u1=1:NG+1;
for loop=1:Nunc
    Big_W0T(u1,v1) = mat_W0.';
    Big_WT(u1,v1) = mat_W.';
    Big_JCds0(v1,u1) = mat_JC_ds0;
    Big_JCds(v1,u1) = mat_JC_ds;
    
    v1=v1+N+1; u1=u1+NG+1;
end 

u_vec = -(Big_W0T*Phi_u0 - Big_WT*Psi_u)\(Big_W0T*Phi_I0 - Big_WT*Psi_I);

%abs(det(Big_W0T*Phi_u0 - Big_WT*Psi_u))

% - calc Phi, Psi, w and wn from this

Phi = Phi_I0 + Phi_u0*u_vec;
Psi = Psi_I + Psi_u*u_vec;

Phi_n = Big_JCds0*u_vec;
Psi_n = Big_JCds*u_vec;

Phi=reshape(Phi,N+1,Nunc);
Psi=reshape(Psi,N+1,Nunc);
Phi_n=reshape(Phi_n,N+1,Nunc);
Psi_n=reshape(Psi_n,N+1,Nunc);

w_coeffs = w_I + w_u*u_vec;
w_n_coeffs = wn_I + wn_u*u_vec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Plot fns u, v, w, w_n, w_s - %%
if 0
figure(fig); hold on
% col='r';
for loop_s=-Nterms:Nterms
    chi_mat(Nterms+1+loop_s,:) = exp(i*loop_s*pi*svec/LC);
end
chi_mat = chi_mat/2/LC;
% for loop_s=-Nterms:Nterms
%     chi_s_mat(Nterms+1+loop_s,:) = (i*loop_s*pi/LC)*exp(i*loop_s*pi*svec/LC);
% end
% chi_s_mat = chi_s_mat/2/LC;
subplot(2,3,1); hold on
title('u')
plot(svec/LC, mat_W0*abs(Phi*chi_mat),col,svec/LC, mat_W*abs(Psi*chi_mat),col);
subplot(2,3,2); hold on
title('v')
plot(svec/LC, inv(mat_JC_ds0)*abs(Phi_n*chi_mat),col,svec/LC, inv(mat_JC_ds)*abs(Psi_n*chi_mat),col);
subplot(2,3,3); hold on
title('w')
plot(svec/LC, abs(w_coeffs.'*chi_mat),col);
subplot(2,3,5); hold on
title('w_n')
plot(svec/LC, abs(w_n_coeffs.'*chi_mat),col);
subplot(2,3,4); hold on
title('w_s')
plot(svec/LC, abs((i*pi/LC)*w_coeffs.'*diag([-Nterms:Nterms])*chi_mat),col);
clear chi_mat col
elseif 0
figure(fig); hold on
for loop_s=-Nterms:Nterms
    chi_mat(Nterms+1+loop_s,:) = exp(i*loop_s*pi*svec/LC);
end
subplot(2,3,1); hold on
title('u')
plot(svec/LC, mat_W*abs(Psi*chi_mat),col);
subplot(2,3,2); hold on
title('w')
plot(svec/LC, abs(w_coeffs.'*chi_mat),col);
subplot(2,3,3); hold on
title('\^{w}')
plot(svec/LC, abs(what_coeffs.'*chi_mat),col);
subplot(2,3,4); hold on
title('u_s')
plot(svec/LC, mat_W*abs((i*pi/LC)*Psi*diag([-Nterms:Nterms])*chi_mat),col);
subplot(2,3,5); hold on
title('w_s')
plot(svec/LC, abs((i*pi/LC)*w_coeffs.'*diag([-Nterms:Nterms])*chi_mat),col);
subplot(2,3,6); hold on
title('\^{w}_s')
plot(svec/LC, abs((i*pi/LC)*what_coeffs.'*diag([-Nterms:Nterms])*chi_mat),col);
clear chi_mat col    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Edge = Psi(1,:);
Edge = sum(abs(Edge/2/LC).^2);
Edge=Edge/2/LC;

%% get the coefficients of del2w & del2w_n:
what_coeffs=( MB1*w_coeffs+MB2*w_n_coeffs );
what_n_coeffs=( MS1*w_coeffs+MS2*w_n_coeffs );

%% CALCULATE R_p vectors:
Rp=zeros(N+1,Nunc);
Sp=Rp;
ddt=zeros(N+1,N+1);
ddt(1)=1;
%%
a_hat=-Phi;
b_hat=Phi_n;

%% calc bessel functions using polar coordinates relative
%% to average of outer contour:
% xy_av=[mean(xvec),mean(yvec)];
% [rr_pol,th_pol]=...
%     GEN_polar_coords( xvec-xy_av(1),yvec-xy_av(2) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% LB(13.11.08)  
 [rr_pol,th_pol]=...
    GEN_polar_coords( xvec,yvec );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
jj=-(1+Nmodes_out):Nmodes_out+1;
jj0=1:(2*Nmodes_out+1);
for lp=1:length(jj)
 besj(:,lp)=besselj(jj(lp),kap0*rr_pol);
end
dr_besj=kap0/2*( besj(:,jj0)-besj(:,jj0+2) );
besj=besj(:,jj0+1);
%%
pp=jj(jj0+1);
Exp=exp( -i*th_pol*pp );
dExp=exp( -i*th_pol*pp )*diag(-i*pp );
%%
Jp=besj.*Exp;
dx_Jp=diag( cos(th_pol) )*( dr_besj.*Exp ) + ...
        diag( -sin(th_pol)./rr_pol )*( besj.*dExp );
dy_Jp=diag( sin(th_pol) )*( dr_besj.*Exp ) + ...
        diag( cos(th_pol)./rr_pol )*( besj.*dExp );
dn_Jp=diag( sin(th_vec) )*dx_Jp + ...
        diag( -cos(th_vec) )*dy_Jp;

%% integrate products with basis functions:
% Jmat0=conj(IP)*Jp;
% Jmat1=conj(IP)*dn_Jp;
%Sp=a_hat*Jmat0+b_hat*Jmat1;
% Rp=1/4i*ddt*Sp;
%%%%%%%%%%%%%%%%%%%%%%%%%
 %% LB (13-14.11.08) and (27.01.08) %%
Jmat0=conj(IP)*Jp;
Jmat1=conj(IP)*dn_Jp;
Sp=a_hat*Jmat1+b_hat*Jmat0;
Rp=-1/4i*ddt*Sp;
% - note the minus sign !
% Rp = Rp(1,:).*exp(pp*pi/2i); % - do this in the main fn !!!!
%%%%%%%%%%%%%%%%%%%%%%%%%

%% - LB (10.11.08):
% 
% a_tilde = a_tilde./(2*LC);
% b_tilde = b_tilde./(2*LC);

Psi_hat = [Psi;w_coeffs.';what_coeffs.'];
Psi_n_hat = [Psi_n;w_n_coeffs.';what_n_coeffs.'];
