function [U_nonsing,dn0U,fac_ds0U,U_nonsing_end]=...
               POOL_hankel_stuff(matrices,vectors,LC,end_points)
%% CALL: [U,dn0U,fac_ds0U]=POOL_hankel_stuff(matrices,vectors)
%%
%% INPUTS:
%% matrices={R,varTH,TH,R_end};
%% vectors={k_in,d2xy_ds2,svec/LC};
%%
%% OUTPUTS:
%% U_nonsing=U-log(R)/2/pi,
%%  U=1/4/pi*H_0(k_in*R);
%% dn0U=\pa_{n_0}U;
%% fac_ds0U=( 1-\cos(\Delta) )*\pa_{s_0}U,
%%  \Delta=\theta(s)-\theta(s_0);

R=matrices{1};
varTH=matrices{2};
TH=matrices{3};
%%
k_in=vectors{1};
N=length(k_in)-1;
d2xy_ds2=vectors{2};
N=length(k_in)-1;
svec_on_LC=vectors{3};%%for testing:
%  [S0oLC,SoLC]=meshgrid(svec_on_LC,svec_on_LC);
%  DSoLC=abs(SoLC-S0oLC);
%%

TH0=TH';
th_vec=TH(:,1);
CD=cos(TH-TH0);
SD=sin(TH-TH0);

CVD=cos(TH-varTH);
SVD=sin(TH-varTH);
CVD0=cos(TH0-varTH);
SVD0=sin(TH0-varTH);

%% NB varTH is bad at R=0,
%% so need to make some corrections to diagonal elements:
CVD=CVD-diag(diag(CVD))+eye(length(TH));
SVD=SVD-diag(diag(SVD));
CVD0=CVD0-diag(diag(CVD0))+eye(length(TH));
SVD0=SVD0-diag(diag(SVD0));

%%
jzero=find(R==0);
jnz=find(R);
logR=0*R;
logR(jnz)=log(R(jnz));
RlogR=R.*logR;
%%
U_nonsing=zeros([size(R),N+1]);
dn0U=U_nonsing;
fac_ds0U=dn0U;
%%
A_in=1/4i+0*k_in;


%% LIMIT OF H_0 - 2i/pi*log(R):
sing_coeff=2i/pi;%% NB sing_coeff/4i=1/2/pi
Hank0_lim=1-sing_coeff*( log(2)+psi(1) );%%check!!!!!!!!

%% LIMIT OF -sin(\Theta_0-\theta_\Delta)/r_\Delta
%% (COMMON TO INNER AND OUTER ROUTINES):
LIM=.5*( d2xy_ds2(:,2).*cos(th_vec) - ...
           d2xy_ds2(:,1).*sin(th_vec) );

%% NOW CALCULATE THE MATRICES:
if nargout>3
  R_end=matrices{4};
  for j=0:N
    %% U's:
    Hank0=besselh( 0,1,k_in(j+1)*R(jnz) ) - sing_coeff*log(R(jnz));
    dummy=0*R+A_in(j+1)*( Hank0_lim+sing_coeff*log(k_in(j+1)) );
    dummy(jnz)=A_in(j+1)*Hank0;
    U_nonsing(:,:,j+1)=dummy;
%    U_in_mn(:,:,j+1)=IP*U*IP'+log_corr_mn;
    %% sing subtracted is 1/2/pi*log(2|s-s_0|),
    %% so have added a correction term.

    U_nonsing_end(:,j+1)=besselh( 0,1,k_in(j+1)*R_end ) -...
          sing_coeff*log( 2*abs(1-svec_on_LC) );

    %% normal derivatives (\pa_{n_0}) of the U's:
    dR_U=-k_in(j+1)*A_in(j+1)*besselh( 1,1,k_in(j+1)*R(jnz) );
    Mat0=A_in(j+1)*sing_coeff*diag(LIM);
    Mat0(jnz)=-SVD0(jnz).*dR_U;
    dn0U(:,:,j+1)=Mat0;

    %% tangential derivatives (\pa_{s_0}) of the U's
    %% ( multiplied by 1-\cos(\Delta) ):
    Mat0=zeros(size(R));
    Mat0(jnz)=-( 1-CD(jnz) ).*CVD0(jnz).*dR_U;
    fac_ds0U(:,:,j+1)=Mat0;
  end
else
  for j=0:N
    %% U's:
    Hank0=besselh( 0,1,k_in(j+1)*R(jnz) ) - sing_coeff*log(R(jnz));
    dummy=0*R+A_in(j+1)*( Hank0_lim+sing_coeff*log(k_in(j+1)) );
    dummy(jnz)=A_in(j+1)*Hank0;
    U_nonsing(:,:,j+1)=dummy;
    %% sing subtracted is 1/2/pi*log(2|s-s_0|),
    %% so have added a correction term.

    %% normal derivatives (\pa_{n_0}) of the U's:
    dR_U=-k_in(j+1)*A_in(j+1)*besselh( 1,1,k_in(j+1)*R(jnz) );
    Mat0=A_in(j+1)*sing_coeff*diag(LIM);
    Mat0(jnz)=-SVD0(jnz).*dR_U;
    dn0U(:,:,j+1)=Mat0;

    %% tangential derivatives (\pa_{s_0}) of the U's
    %% ( multiplied by 1-\cos(\Delta) ):
    Mat0=zeros(size(R));
    Mat0(jnz)=-( 1-CD(jnz) ).*CVD0(jnz).*dR_U;
    fac_ds0U(:,:,j+1)=Mat0;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0%% CHECK THE INNER HANKEL FUNCTIONS' EXPANSION:
  disp('checking inner H_0s:');
  jchq=0;
  j=jchq;
  Nint=size(R,1);
  tstA=[A_in(j+1),1/4i]
  %%
  U_j=A_in(j+1)*besselh(0,1,k_in(j+1)*R);
  U_j_mn=U_in_mn(:,:,j+1);

  nchk=11%GEN_rand(Nint)
  Nterms=(size(IP,1)-1)/2;
  nnfm=(-Nterms:Nterms)';
  Am=U_j_mn*exp(-i*nnfm*pi*svec_on_LC(nchk));
  plot( svec_on_LC,U_j(:,nchk) ), hold on;
  Ujap=GEN_interp_exp(svec_on_LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot( svec_on_LC,Ujap,'--r' ), hold off;

  if 1%% IF CIRCULAR GEOMETRY:
    radius=10;
    tst_circ=[diag(U_j_mn),...
                A_in(j+1)*besselj(nnfm,k_in(j+1)*radius).*...
                   besselh(nnfm,1,k_in(j+1)*radius)]
  end

elseif 0%% CHECK THE NORMAL DERIVATIVE
        %% OF THE INNER HANKEL FUNCTIONS' EXPANSION:
  disp('checking norm derivs of inner H_0s:');
  jchq=0;
  j=jchq;
  Nint=size(R,1);
  %tstA=[A_in(j+1),1/4i]
  %%
  jz=find(R==0);
  dn0U_j=0*R;
  dn0U_j(jnz)=k_in(j+1)*A_in(j+1)*...
     besselh(1,1,k_in(j+1)*R(jnz)).*SVD0(jnz);
  dn0U_j=dn0U_j+A_in(j+1)*sing_coeff*diag(LIM);
  dn0U_j_mn=dn0U_in_mn(:,:,j+1);

  nchk=GEN_rand(Nint)
  Nterms=(size(IP,1)-1)/2;
  nnfm=(-Nterms:Nterms)';
  Am=dn0U_j_mn*exp(-i*nnfm*pi*svec_on_LC(nchk));
  plot( svec_on_LC,dn0U_j(:,nchk) ), hold on;
  Ujap=GEN_interp_exp(svec_on_LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot( svec_on_LC,Ujap,'--r' ), hold off;

  if 1%% IF CIRCULAR GEOMETRY:
    %k_in(j+1)=2;
    radius=10;
    dJ=-.5*( besselj(nnfm+1,k_in(j+1)*radius)-...
                       besselj(nnfm-1,k_in(j+1)*radius) );
    dH=-.5*( besselh(nnfm+1,1,k_in(j+1)*radius)-...
                       besselh(nnfm-1,1,k_in(j+1)*radius) );
    wronsk=2i/pi/radius;
    CC=A_in(j+1)*.5*wronsk;%[CC,1/4/pi/radius]
    tst_circ=[diag(dn0U_j_mn), diag( IP*dn0U_j*IP' )...
                 CC+A_in(j+1)*k_in(j+1)*dJ.*...
                      besselh(nnfm,1,k_in(j+1)*radius),...
                 -CC+A_in(j+1)*k_in(j+1)*dH.*...
                      besselj(nnfm,k_in(j+1)*radius)]
  else
    [ diag(dn0U_j_mn), diag( IP*dn0U_j*IP' ) ]
  end

elseif 0%% CHECK THE OUTER HANKEL FUNCTIONS' EXPANSION:
  disp('checking outer H_0s:');
  jchq=0;
  j=jchq;
  Nint=size(R,1);
  tstA=[A_out(j+1),1/4i]
  %%
  U_j=A_out(j+1)*besselh(0,1,k_out(j+1)*R);
  U_j_mn=U_out_mn(:,:,j+1);

  nchk=GEN_rand(Nint)
  Nterms=(size(IP,1)-1)/2;
  nnfm=(-Nterms:Nterms)';
  Am=U_j_mn*exp(-i*nnfm*pi*svec_on_LC(nchk));
  plot( svec_on_LC,U_j(:,nchk) ), hold on;
  Ujap=GEN_interp_exp(svec_on_LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot( svec_on_LC,Ujap,'--r' ), hold off;

  if 1%% IF CIRCULAR GEOMETRY:
    radius=10;
    tst_circ=[diag(U_j_mn),...
                A_out(j+1)*besselj(nnfm,k_out(j+1)*radius).*...
                   besselh(nnfm,1,k_out(j+1)*radius)]
  end

elseif 0%% CHECK THE NORMAL DERIVATIVE
        %% OF THE OUTER HANKEL FUNCTIONS' EXPANSION:
  disp('checking norm derivs of outer H_0s:');
  jchq=0;
  j=jchq;
  Nint=size(R,1);
  %tstA=[A_out(j+1),1/4i]
  %%
  jz=find(R==0);
  dn0U_j=0*R;
  dn0U_j(jnz)=k_out(j+1)*A_out(j+1)*...
     besselh(1,1,k_out(j+1)*R(jnz)).*SVD0(jnz);
  dn0U_j=dn0U_j+A_out(j+1)*sing_coeff*diag(LIM);
%tst_dnU=dn0U_j(1:14,1:14)
  dn0U_j_mn=dn0U_out_mn(:,:,j+1);

  nchk=GEN_rand(Nint)
  Nterms=(size(IP,1)-1)/2;
  nnfm=(-Nterms:Nterms)';
  Am=dn0U_j_mn*exp(-i*nnfm*pi*svec_on_LC(nchk));
  plot( svec_on_LC,dn0U_j(:,nchk) ), hold on;
  Ujap=GEN_interp_exp(svec_on_LC,Am);
  %[Lap,log_corr(:,nchk)]
  plot( svec_on_LC,Ujap,'--r' ), hold off;

  if 1%% IF CIRCULAR GEOMETRY:
    %k_out(j+1)=2;
    radius=10;
    dJ=-.5*( besselj(nnfm+1,k_out(j+1)*radius)-...
                       besselj(nnfm-1,k_out(j+1)*radius) );
    dH=-.5*( besselh(nnfm+1,1,k_out(j+1)*radius)-...
                       besselh(nnfm-1,1,k_out(j+1)*radius) );
    wronsk=2i/pi/radius;
    CC=A_out(j+1)*.5*wronsk;%[CC,1/4/pi/radius]
    tst_circ=[diag(dn0U_j_mn), diag( IP*dn0U_j*IP' )...
                 CC+A_out(j+1)*k_out(j+1)*dJ.*...
                      besselh(nnfm,1,k_out(j+1)*radius),...
                 -CC+A_out(j+1)*k_out(j+1)*dH.*...
                      besselj(nnfm,k_out(j+1)*radius)]
  elseif 0
    [ diag(dn0U_j_mn), diag( IP*dn0U_j*IP' ) ]
  end
end