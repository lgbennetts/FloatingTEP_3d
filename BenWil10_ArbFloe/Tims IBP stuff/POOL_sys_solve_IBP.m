function [a_tilde,b_tilde,Rp,LC,svec, dtheta_ds]=... %XTRA_OUTPUTS]=...
             POOL_sys_solve_IBP(matrices,vectors,...
               parameters,edge_vars,N,Nterms,Nmodes_out, Nint)
%% CALL: POOL_sys_solve_smooth(matrices,vectors,parameters,...
%%        edge_vars,Nterms)
%% matrices={Pminus,Pplus,CMAT,Aplus}
%% vectors={k_in,k_out,dAminus}
%% parameters={alf,beta,sigma,nu,theta_inc in degrees}
%% edge_vars={fxn_handle,fxn_parameters,srt}
%% srt is best seen by looking below in the test inputs

USE_EXP=1;%% use exponentials for basis,
          %% else use cos's & sin's.
NORMAL_POINTING_OUT=1;
NORM_FAC=(2*NORMAL_POINTING_OUT-1);

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
%%
%Nint=300;
if 0
  [tt,ww]=OP_numint_legendre(Nint);%% use Gauss-legendre points
                          %% and weights as everything smooth.
  [xyvecs,ds_dt,th_vec,dtheta_ds,d2s_dt2,d2theta_ds2,...
     d2xy_ds2,LC,svec]=CRK_get_rsdtheta(edge_vars,tt);
  ww=ww.*ds_dt/LC;
  %% (1/LC)\int_{-LC}^{LC}ds=(1/LC)\int_{-1}^{1}(ds/dt)dt
  Npts=Nint;
elseif 1
  %% use 'constant panel' integration scheme;
  %% no of points is 2*Nint;
  %% NB this int scheme exactly integrates
  %%  Fourier modes with |index|<=(Nint-1)
  %%  => need (Nint-1)>=Nterms;
  %% ALSO NB integrates with respect to s, instead of the
  %%  curve parameter t.
  Nint=max(2*75,Nterms+1);
  Npts=2*Nint
  [svec_on_LC,ww]=GEN_numint_exp(Nint);
  %% (1/LC)\int_{-LC}^{LC}ds=\int_{-1}^{1}d(s/LC)
  [xyvecs,ds_dt,th_vec, dtheta_ds,...
    d2s_dt2,d2theta_ds2,d2xy_ds2,LC,t_quad]=...
      CRK_get_rsdtheta_NR( edge_vars,svec_on_LC );
  svec=svec_on_LC*LC;
  save('GEOMSTUFF.mat','xyvecs','ds_dt','th_vec',...
        'dtheta_ds','d2s_dt2','d2theta_ds2','d2xy_ds2',...
         'LC','t_quad','svec','ww');
else
  load GEOMSTUFF;
  Npts=length(svec)
  Nint=Npts/2;
%  [svec(1:5)/LC t_quad(1:5);...
%    svec(end-4:end)/LC t_quad(end-4:end)]
%  dx=xyvecs(:,1)-xyvecs(:,end);
%  dis=sqrt(dx'*dx)
end

xvec=xyvecs(1,:).';
yvec=xyvecs(2,:).';
[X0,X]=meshgrid(xvec,xvec);
[Y0,Y]=meshgrid(yvec,yvec);
[R,varTH]=GEN_polar_coords(X-X0,Y-Y0);
[TH0,TH]=meshgrid(th_vec,th_vec);
jnz=find(R);
geom_stuff={xyvecs,ds_dt,th_vec,dtheta_ds,d2s_dt2,d2theta_ds2,...
             d2xy_ds2,t_quad,svec/LC};

if 0%% CHECK SHAPE & TANGENT/NORMAL VECTORS:
  %dtheta_ds,1/radius
  plot(xvec,yvec), hold on;
  svectors=[cos(th_vec),sin(th_vec)].';
  nvectors=[sin(th_vec),-cos(th_vec)].';
  dj=5;
  for j=1:dj:Npts
    svector_x=[0;svectors(1,j)];
    svector_y=[0;svectors(2,j)];
    plot(xvec(j)+svector_x,yvec(j)+svector_y,'m');
    %%
    nvector_x=[0;nvectors(1,j)];
    nvector_y=[0;nvectors(2,j)];
    plot(xvec(j)+nvector_x,yvec(j)+nvector_y,'g');
  end
  daspect([1 1 1]);
  hold off;
  return;
end

%% GET INNER PRODUCT MATRICES:
%ds_dt,pause
[IP,hn,Exp]=GEN_inprod_exp(svec/LC,ww,Nterms); %_vOld
nn_unc=(-Nterms:Nterms).';
kunc=pi/LC*nn_unc;
Ddiff=diag(i*kunc);
ip_stuff={IP,hn,Exp};
%  figure,plot(svec/LC,[dtheta_ds,...
%    GEN_interp_exp(svec/LC,XTRA_OUTPUTS{2})]); return
%hn=[2*LC;LC+0*(1:Nunc-1)'];

k_wtr=vectors{1};
k_ice=vectors{2};
N=length(k_wtr)-1;
%  matrices4U={R,varTH,TH};
%  vectors4U={k_in,k_out,d2xy_ds2,svec/LC};

%% GET LOG CORRECTION:
%% need to add 1/2\pi*log(R) singularity into U's
%% ( this is the same for all the k's ):
[S0,SS]=meshgrid(svec,svec);
Exp_term=abs( 1-exp(pi*i/LC*(SS-S0)) );

%% coefficients for 1/(2*pi)*log[1-exp(i*(s-s_0)/LC)]:
ln_coeff=-1/4/pi./(1:Nterms).';
ln_coeff=[flipud(ln_coeff);0;ln_coeff];
log_corr=diag( 0*ds_dt - 1/2/pi*log(pi/LC) );
log_corr(jnz)=1/2/pi*log( R(jnz)./Exp_term(jnz) );
log_corr_mn=diag(ln_coeff)+IP*log_corr*IP';
log_stuff={log_corr_mn};
%%
%% GET MATRICES FOR INNER PROBLEM:
Pminus=matrices{1};
PminusT=Pminus.';
PminusT_inv=inv(PminusT);
dAminus=vectors{3};
%%
M00left=PminusT;
M00right=PminusT_inv;
M01left=M00left;
M01right=-diag(1./dAminus)*Pminus;
Mmats_u={M00left, M00right,M01left,M01right};
%%
[KMAT0,MDIAG0] = POOL_kernelmat_wtr_diag_exp(k_wtr,Mmats_u,...
                   geom_stuff,LC,ip_stuff,log_stuff);
KMAT_in=[ MDIAG0-NORM_FAC*KMAT0,zeros( (N+1)*Nunc,2*Nunc ) ];
%% NB this matrix has to be padded out with zeros to allow for
%%  the (w,w_n) columns;
%% ALSO NB the minus sign is to ensure we solve
%%  u/2 - L[u,v]=f,
%%   since the POOL_kernelmat_wtr_diag_exp.m matrix
%%    corresponds to L[u,v]-u/2;
FORCING_in=KMAT_in(:,end);

%% GET MATRICES FOR OUTER PROBLEM:
Pplus=matrices{2};
Pplus_tilde=eye(N+3);
Pplus_tilde(1:N+1,1:N+1)=Pplus;
PplusT_tilde=Pplus_tilde.';
invPplusT_tilde=inv(PplusT_tilde);
%%
CMAT=matrices{3};
invCMAT=inv(CMAT);
%%
Aplus=matrices{4};
invAplus_tilde=eye(N+3);
invAplus_tilde(1:N+1,1:N+1)=inv(Aplus);
%%
M11left=PplusT_tilde*CMAT;
M11right=-invCMAT*invAplus_tilde*Pplus_tilde;
M10left=M11left;
M10right=inv(M11left);
Mmats_u={M10left,M10right,M11left,M11right};
%%
%% BENDING MOMENT CONDITION:
%% del2w_coeffs=MB1*w_coeffs +MB2*dnw_coeffs 
MB1=diag(-bet0*kunc.^2);%% from \beta_0*\pa_s^2w
MB2=diag(hn*LC)*IP*diag(bet0*dtheta_ds)*Exp*diag( (1/LC)./hn );
%% from bet0*dtheta_ds*\pa_nw

%% SHEARING CONDITION:
%% dn_del2w_coeffs=MS1*w_coeffs +MS2*dnw_coeffs
MS1=Ddiff*MB2*Ddiff;%% from \beta_0\pa_s*dtheta_ds*\pa_sw
MS2=-MB1;%% from -\beta_0*\pa_s^2\pa_nw

if 0
  edge_mats={MB1,MB2,MS1,MS2};
  [KMAT1,MDIAG1]=POOL_kernelmat_ice_diag_exp_v0(k_ice,...
     Mmats_u,geom_stuff,LC,ip_stuff,log_stuff,edge_mats);
else
  Mwns1=diag(bet0*hn*LC)*IP*diag(dtheta_ds)*Exp*...
    diag( (i*kunc/LC)./hn );% MS1, Ddiff*Mwns1
  Mwns2=-bet0*Ddiff;% MS2, Ddiff*Mwns2
  edge_mats={MB1,MB2,Mwns1,Mwns2};
  [KMAT1,MDIAG1]=POOL_kernelmat_ice_diag_exp(k_ice,...
     Mmats_u,geom_stuff,LC,ip_stuff,log_stuff,edge_mats);
end
KMAT_out=MDIAG1+NORM_FAC*KMAT1;

FORCING_out=0*KMAT_out(:,1);
kap0=k_ice(1);
ff=exp(i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
fm=IP*ff;
for j=0:N+2
  jvec0=(1:Nunc)+j*Nunc;
  FORCING_out(jvec0)=M11left(j+1,1)*fm;
end

%% SOLVE INTEGRAL EQUATION AND REARRANGE OUTPUTS:
FORCING=[FORCING_in;FORCING_out];
KMAT=[KMAT_in;KMAT_out];
ab_tilde=KMAT\FORCING;
%%
a_tilde=zeros(N+3,Nunc);
b_tilde=zeros(N+3,Nunc);
for j=0:N+1
  jv_a=(1:Nunc)+2*j*Nunc;
  jv_b=jv_a+Nunc;
  a_tilde(j+1,:)=ab_tilde(jv_a).';
  b_tilde(j+1,:)=ab_tilde(jv_b).';
end

%% GET THE COEFFICIENTS OF del2w & del2w_n:
w_coeffs=a_tilde(N+2,:).';
w_n_coeffs=b_tilde(N+2,:).';
a_tilde(N+3,:)=( MB1*w_coeffs+MB2*w_n_coeffs ).';
b_tilde(N+3,:)=( MS1*w_coeffs+MS2*w_n_coeffs ).';

%% CALCULATE R_p vectors:
Rp=zeros(N+3,Nunc);
Sp=Rp;
ddt=zeros(N+3,N+3);
ddt(1)=1;
%%
a_hat=( -NORM_FAC*M10right*a_tilde );
b_hat=( -NORM_FAC*M11right*b_tilde );
%XTRA_OUTPUTS={svec/LC,MB2,hn};

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
besj=besselj(jj,kap0*rr_pol);
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
Rp=-1/4i*CMAT*ddt*Sp;
% - note the minus sign !
% Rp = Rp(1,:).*exp(pp*pi/2i); % - do this in the main fn !!!!
%%%%%%%%%%%%%%%%%%%%%%%%%

%% - LB (10.11.08):
% 
% a_tilde = a_tilde./(2*LC);
% b_tilde = b_tilde./(2*LC);