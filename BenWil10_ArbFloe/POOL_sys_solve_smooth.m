function [a_tilde,b_tilde,Rp,LC]=...
             POOL_sys_solve_smooth(matrices,vectors,...
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
xvec=xyvecs(1,:).';
yvec=xyvecs(2,:).';
[X0,X]=meshgrid(xvec,xvec);
[Y0,Y]=meshgrid(yvec,yvec);
[R,varTH]=GEN_polar_coords(X-X0,Y-Y0);
[TH0,TH]=meshgrid(th_vec,th_vec);
jnz=find(R);

if 0%% CHECK SHAPE & TANGENT/NORMAL VECTORS:
  %dtheta_ds,1/radius
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
%[U_in_mn,U_out_mn,dn0U_in_mn,dn0U_out_mn]=...
%  POOL_hankel_derivs(matrices4U,vectors4U,ipstuff);
%U_in_mn(:,:,1)
%dn0U_out_mn(:,:,3)
[U_in_mn,U_out_mn,dn0U_in_mn,dn0U_out_mn]=...
POOL_hankel_derivs2(matrices4U,vectors4U,ipstuff);
%dn0U_out_mn(:,:,3)

%% CALC KERNEL MATRIX FROM INNER PROBLEM:
Pminus=matrices{1};
PminusT=Pminus.';
PminusT_inv=inv(PminusT);
dAminus=vectors{3};
%%
M00left=PminusT;
M00right=PminusT_inv;
%%
M01left=M00left;
M01right=-diag(1./dAminus)*Pminus;
%%
Jvec=(1:Nunc)';
jv_w=Jvec+(N+1)*Nunc;
jv_del2w=jv_w+Nunc;
jh_w=Jvec+2*(N+1)*Nunc;
jh_dnw=jh_w+Nunc;
%%
KMAT_in=zeros((N+1)*Nunc,(2*N+4)*Nunc);
forcing_in=KMAT_in(:,1);%% this is zero on the inside,
                        %% so leave as is.
Mdiag=diag(.5./hn);%% NB factor 1/2;
   %% this matrix goes on the diagonal of the main matrix
   %% - hn is because we take norm of chi_n so fac of 2LC (LB 11.11.08)

for j=0:N
  jvec0=Jvec+j*Nunc;
  jvec1=Jvec+2*j*Nunc;
  KMAT_in(jvec0,jvec1) = KMAT_in(jvec0,jvec1) + Mdiag;
  for r=0:N
    jvec1=Jvec+2*r*Nunc;
    jvec2=jvec1+Nunc;
    for p=0:N
      KMAT_in(jvec0,jvec1) = KMAT_in(jvec0,jvec1) - ...
        M00left(j+1,p+1)*M00right(p+1,r+1)*dn0U_in_mn(:,:,p+1);
      KMAT_in(jvec0,jvec2) = KMAT_in(jvec0,jvec2) - ...
        M01left(j+1,p+1)*M01right(p+1,r+1)*U_in_mn(:,:,p+1);
    end
  end
end

%% CALC KERNEL MATRIX FROM OUTER PROBLEM:
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
M11right=invCMAT*invAplus_tilde*Pplus_tilde;
%%
M10left=M11left;
M10right=-inv(M11left);
%%
KMAT_out=zeros((N+3)*Nunc,(2*N+4)*Nunc);

%% FORCING TERMS:
forcing_out=KMAT_out(:,1);
kap0=k_out(1);
ff=exp(i*kap0*(xvec*cos(theta_inc)+yvec*sin(theta_inc)));
%[ff,exp(i*kap0*xvec)],pause
%fm=IP*ff;
%%%%%%% LB (12.11.08) %%%%%%%%
% fm=conj(IP)*ff;
fm=IP*ff;
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
 %% - This part does the lhs bending moment manip (LB 11.11.08)  
KMAT_out(jv_del2w,jh_w)=Mdiag*MB1;
KMAT_out(jv_del2w,jh_dnw)=Mdiag*MB2;

%% matrices to apply shearing condition:
%% dn_del2w_coeffs=MS1*w_coeffs +MS2*dnw_coeffs
MS1=Ddiff*MB2*Ddiff;%% from \beta_0\pa_s*dtheta_ds*\pa_sw
%  %%%%%%% LB (12.11.08) %%%%%%%%
%  MS1=-Ddiff*MB2*Ddiff;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
MS2=-MB1;%% from -\beta_0*\pa_s^2\pa_nw

%[Nunc,size(KMAT_out),2*(N+2)*Nunc]
for j=0:N+2
  jvec0=Jvec+j*Nunc;
  forcing_out(jvec0)=M11left(j+1,1)*fm;
  %%
  jvec1=Jvec+2*j*Nunc;%j,jvec1'
  if j<N+2
    KMAT_out(jvec0,jvec1) = KMAT_out(jvec0,jvec1) + ...
      diag(1/2./hn);%% NB factor 1/2,
      %% also NB need to do something different for j=N+2
%  KMAT_out(1:4*Nunc,1:4*Nunc),pause
%    else
%      KMAT_out(1:4*Nunc,1:4*Nunc)
  end
  for r=0:N+1
    jvec1=Jvec+2*r*Nunc;
    jvec2=jvec1+Nunc;
    for p=0:N+2
      KMAT_out(jvec0,jvec1) = KMAT_out(jvec0,jvec1) - ...
        M10left(j+1,p+1)*M10right(p+1,r+1)*...
          dn0U_out_mn(:,:,p+1);
      KMAT_out(jvec0,jvec2) = KMAT_out(jvec0,jvec2) - ...
        M11left(j+1,p+1)*M11right(p+1,r+1)*U_out_mn(:,:,p+1);
    end
  end
  for r=N+2
    for p=0:N+2
      KMAT_out(jvec0,jh_w) = KMAT_out(jvec0,jh_w) - ...
        M10left(j+1,p+1)*M10right(p+1,r+1)*...
           dn0U_out_mn(:,:,p+1)*MB1;
      KMAT_out(jvec0,jh_dnw) = KMAT_out(jvec0,jh_dnw) - ...
        M10left(j+1,p+1)*M10right(p+1,r+1)*...
           dn0U_out_mn(:,:,p+1)*MB2;
      %%
      KMAT_out(jvec0,jh_w) = KMAT_out(jvec0,jh_w) - ...
        M11left(j+1,p+1)*M11right(p+1,r+1)*...
           U_out_mn(:,:,p+1)*MS1;
      KMAT_out(jvec0,jh_dnw) = KMAT_out(jvec0,jh_dnw) - ...
        M11left(j+1,p+1)*M11right(p+1,r+1)*...
           U_out_mn(:,:,p+1)*MS2;
    end
  end
end

%% SOLVE INTEGRAL EQUATION AND REARRANGE OUTPUTS:
FORCING=[forcing_in;forcing_out];
KMAT=[KMAT_in;KMAT_out];
% ab_tilde=-FORCING\KMAT;

 %% - LB (11.11.08)
 
ab_tilde=KMAT\FORCING; 

%%
a_tilde=zeros(N+3,Nunc);
b_tilde=zeros(N+3,Nunc);
for j=0:N+1
  jv_a=Jvec+2*j*Nunc;
  jv_b=jv_a+Nunc;
  a_tilde(j+1,:)=ab_tilde(jv_a).';
  b_tilde(j+1,:)=ab_tilde(jv_b).';
end

%% get the coefficients of del2w & del2w_n:
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
a_hat=( M10right*a_tilde );
b_hat=( M11right*b_tilde );

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