function [U_in_mn,U_out_mn,dn0U_in_mn,dn0U_out_mn]=...
   POOL_hankel_derivs2(matrices,vectors,ipstuff)
%%matrices={R,varTH,TH};
%%vectors={k_in,k_out,d2xy_ds2,svec/LC}

U_in_mn = [];U_out_mn=[]; dn0U_in_mn = []; dn0U_out_mn = [];

R=matrices{1};
varTH=matrices{2};
TH=matrices{3};
%%
k_in=vectors{1};
k_out=vectors{2};
d2xy_ds2=vectors{3};
N=length(k_in)-1;
svec_on_LC=vectors{4};%%for testing:
%%
IP=ipstuff{1};
log_corr_mn=ipstuff{2};
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
jnz=find(R);
logR=0*R;
logR(jnz)=log(R(jnz));
RlogR=R.*logR;
%%
U_in=zeros([size(R),N+1]);
dnU_in=U_in;
U_out=zeros([size(R),N+3]);
dnU_out=U_out;
dsU_out=U_out;
%%
A_in=1/4i+0*k_in;
A_out=1/4i+0*k_out;

%% LIMIT OF H_0 - 2i/pi*log(R):
sing_coeff=2i/pi;%% NB sing_coeff/4i=1/2/pi
Hank0_lim=1-sing_coeff*(log(2)+psi(1));
dn0U=zeros(size(R));

%% LIMIT OF -sin(\Theta_0-\theta_\Delta)/r_\Delta
%% (COMMON TO INNER AND OUTER ROUTINES):
LIM=.5*( d2xy_ds2(:,2).*cos(th_vec) - ...
             d2xy_ds2(:,1).*sin(th_vec) );

%% INNER FUNCTIONS:
for j=0:N
  %% U's:
  Hank0=besselh( 0,1,k_in(j+1)*R(jnz) )-sing_coeff*logR(jnz);
  U=0*R+A_in(j+1)*( Hank0_lim+sing_coeff*log(k_in(j+1)) );
  U(jnz)=A_in(j+1)*Hank0;
  U_in_mn(:,:,j+1)=IP*U*IP'+log_corr_mn;
  %% sing subtracted is 1/2/pi*log(R),
  %% so have added a correction term.

  %% normal derivatives of the U's:
  dn0U_j=A_in(j+1)*sing_coeff*diag(LIM);
  dn0U_j(jnz)=k_in(j+1)*A_in(j+1)*...
     besselh(1,1,k_in(j+1)*R(jnz)).*SVD0(jnz);
  dn0U_in_mn(:,:,j+1) = IP*dn0U_j*IP';
end


%% OUTER FUNCTIONS:
for j=0:N+2
  %% U's:
  Hank0=besselh( 0,1,k_out(j+1)*R(jnz) )-sing_coeff*logR(jnz);
  U=zeros(size(R))+A_out(j+1)*(Hank0_lim+sing_coeff*log(k_out(j+1)));
  U(jnz)=A_out(j+1)*Hank0;
  U_out_mn(:,:,j+1)=IP*U*IP'+log_corr_mn;
  %% sing subtracted is 1/2/pi*log(R),
  %% so have added a correction term.

  %% normal derivatives of the U's:
  dn0U_j=A_out(j+1)*sing_coeff*diag(LIM);
  dn0U_j(jnz)=k_out(j+1)*A_out(j+1)*...
     besselh(1,1,k_out(j+1)*R(jnz)).*SVD0(jnz);
%tstdn0U=dn0U_j(1:14,1:14)
  dn0U_out_mn(:,:,j+1) = IP*dn0U_j*IP';
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

