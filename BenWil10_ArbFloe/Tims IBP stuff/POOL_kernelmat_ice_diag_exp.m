function [KMAT,MDIAG]=...
  POOL_kernelmat_ice_diag_exp(k_ice,Mmats_u,...
    geom_stuff,LC,ip_stuff,log_stuff,edge_mats,Mmats_v)
%    POOL_kernelmat_wtr_diag_v0(k_wtr,Mmats_u,Mmats_v,...
%      geom_stuff,LC,end_points,ip_stuff,log_stuff,xtra_ints,Nunc)

%% SMOOTH LEAD ONLY (Fourier modes);
%%
%% OPTIONAL CORRECTION FOR IRREGULAR FREQUENCIES:
%% IF 'Mmats_v' is inputted
%% THEN just solve
%%  (1/2)\bfu=\oint[ \hat{\bfG}_n0\bfu - \hat{\bfG}\bfv ]ds0,
%%   \hat{\bfG}_n0=M00_l*\bfU_n0*M00_r,
%%    -\hat{\bfG}=M01_l*\bfU*M01_r,
%%     where [\bfU]_{jr}=\delta_{jr}/(4i)*H_0(k_jR),
%%      R=|\bfx-\bfx_0|;
%% ELSE also use
%%  (1/2)\bfv=\oint[ \breve{\bfG}_{n,n0}\bfu - \breve{\bfG}_n\bfv ]ds,
%%   \breve{\bfG}_{n,n0}=M10_l*\bfU_{n,n0}*M10_r,
%%    -\breve{\bfG}_n=M11_l*\bfU_n*M11_r,

%% INPUTS:
%% Mmats_u={M00_l,M00_r,M01_l,M01_r};
%% Mmats_v={M10_l,M10_r,M11_l,M11_r};
%% geom_stuff={xyvecs,ds_dt,th_vec,dtheta_ds,d2s_dt2,d2theta_ds2,...
%%              d2xy_ds2,t_quad,svec_on_LC};
%% LC=half-length of curve;
%% ip_stuff={IP,hn,Exp}, where Fourier coefficients
%%  f_m=[ IP*f(x_n) ]_m,
%%   hn=1/2,
%%    [ Exp ]_{mn}=exp(i*k_n*s_m),
%%     k_n=n*pi/LC, n=-Nterms..Nterms;
%% log_stuff={log_corr_mn}
%%  (m,n)-th element of matrix 'log_corr_mn'
%%   is < psi_m | (1/2/pi)*log(R) | psi_n >
%% edge_mats={MB1,MB2,Mwns1,Mwns2} apply the edge conditions by
%%  eliminating \hat{w} & \hat{w_n}:
%%   [\hat{a}_n]=MB1*[a_n]+MB2*[b_n],
%%   [\hat{b}_n]=Mwns1*[a_n]+Mwns2*[\hat{b}_n]

xyvecs=geom_stuff{1};
xvec=xyvecs(1,:).';
yvec=xyvecs(2,:).';
th_vec=geom_stuff{3};
d2xy_ds2=geom_stuff{7};
svec_on_LC=geom_stuff{9};
%%
[X0,X]=meshgrid(xvec,xvec);
[Y0,Y]=meshgrid(yvec,yvec);
[R,varTH]=GEN_polar_coords(X-X0,Y-Y0);
[TH0,TH]=meshgrid(th_vec,th_vec);
SD=sin(TH-TH0);
jnz=find(R);
%%
%  v_end=end_points(:,2);
%  R_end=GEN_polar_coords( xvec-v_end(1),yvec-v_end(2) );

%  USE_XTRA=~isempty(xtra_ints);
%  if USE_XTRA
%    ip0=LC*xtra_ints{1};%% non-smooth lead
%  else
%    ip0=LC*zeros(1,length(xvec));%% smooth lead
%  end

%% INNER PRODUCT MATRICES:
IP=ip_stuff{1};
hn=ip_stuff{2};
Exp=ip_stuff{3};
Nunc=length(hn);
N=length(k_ice)-3;
Nunc_total=(N+2)*(2*Nunc);%% (N+1)*(2*Nunc) [u,v] + 2*Nunc [w,w_n]
Nterms=round((Nunc-1)/2);
kn=pi/LC*(-Nterms:Nterms)';
Ddiff=diag(i*kn);
%%
%matrices={R,varTH,TH,R_end};
matrices={R,varTH,TH};
vectors={k_ice,d2xy_ds2,svec_on_LC};
[U_nonsing,dn0U,fac_ds0U]=...
  POOL_hankel_stuff(matrices,vectors,LC);
log_corr_mn=log_stuff{1};

for j=0:N+2
  U_mn(:,:,j+1)=IP*U_nonsing(:,:,j+1)*IP'+log_corr_mn;
  dn0U_mn(:,:,j+1)=IP*dn0U(:,:,j+1)*IP';
end

%% MATRICES FOR EDGE CONDITIONS:
MB1=edge_mats{1};
MB2=edge_mats{2};
Mwns1=edge_mats{3};
Mwns2=edge_mats{4};

IF_IRREG=(nargin>7);
if ~IF_IRREG
  %% NO CORRECTION FOR IRREGULAR FREQ'S:
  Nrows=(N+3)*Nunc;
  KMAT=zeros(Nrows,Nunc_total);
  MDIAG=KMAT;
  Mdiag=diag( (.5/LC)./hn );
  %%
  M00left=Mmats_u{1};
  M00right=Mmats_u{2};
  M01left=Mmats_u{3};
  M01right=Mmats_u{4};
  %%
  jw_col=(1:Nunc)+(N+1)*(2*Nunc);
  jdnw_col=jw_col+Nunc;
  jwhat_row=(1:Nunc)+(N+2)*Nunc;

  %% THIS ELIMINATES \hat{w} from LHS:
  MDIAG(jwhat_row,jw_col)=Mdiag*MB1;
  MDIAG(jwhat_row,jdnw_col)=Mdiag*MB2;
  %%
  for j=0:N+2
    jvec0=(1:Nunc)+j*Nunc;
    jvec1=(1:Nunc)+2*j*Nunc;%j,jvec1'
    if j<N+2
      MDIAG(jvec0,jvec1) = Mdiag;
      %% NB factor 1/2 in Mdiag;
      %% ALSO NB need to do something different for j=N+2
    end

    %% do (u,v) & (w,w_n) columns:
    for r=0:N+1
      jvec1=(1:Nunc)+2*r*Nunc;%% a^{(r)}_n col
      jvec2=jvec1+Nunc;%% b^{(r)}_n col
      for p=0:N+2
        KMAT(jvec0,jvec1) = KMAT(jvec0,jvec1) + ...
          M00left(j+1,p+1)*M00right(p+1,r+1)*...
            dn0U_mn(:,:,p+1);
        KMAT(jvec0,jvec2) = KMAT(jvec0,jvec2) + ...
          M01left(j+1,p+1)*M01right(p+1,r+1)*U_mn(:,:,p+1);
      end
    end

    %% put (\hat{w},\hat{w}_n) columns into (w,w_n) columns
    %% (USE EDGE CONDITIONS):
    for r=N+2
      for p=0:N+2
        %% deal with \hat{w}:
        KMAT(jvec0,jw_col) = KMAT(jvec0,jw_col) + ...
          M00left(j+1,p+1)*M00right(p+1,r+1)*...
             dn0U_mn(:,:,p+1)*MB1;
        KMAT(jvec0,jdnw_col) = KMAT(jvec0,jdnw_col) + ...
          M00left(j+1,p+1)*M00right(p+1,r+1)*...
             dn0U_mn(:,:,p+1)*MB2;

        %% now deal with \hat{w_n}:
        MUp=Ddiff*U_mn(:,:,p+1)+...
              -IP*fac_ds0U(:,:,p+1).'*IP'+...
                IP*(SD.*dn0U(:,:,p+1).')*IP';
        %% NB these are \pa_sU, [\cos(\th-\th_0)-1]*\pa_sU,
        %%                \sin(\th-\th_0)*\pa_nU

        KMAT(jvec0,jw_col) = KMAT(jvec0,jw_col) + ...
          M01left(j+1,p+1)*M01right(p+1,r+1)*...
             MUp*Mwns1;
        KMAT(jvec0,jdnw_col) = KMAT(jvec0,jdnw_col) + ...
          M01left(j+1,p+1)*M01right(p+1,r+1)*...
             MUp*Mwns2;
      end
    end
  end
else
  %% INCLUDE CORRECTION FOR IRREGULAR FREQ'S:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEED TO FINISH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  KMAT=zeros(Nunc_total,Nunc_total);
  %%
  M00left=Mmats_u{1};
  M00right=Mmats_u{2};
  M01left=Mmats_u{3};
  M01right=Mmats_u{4};
  %%
  M10left=Mmats_v{1};
  M10right=Mmats_v{2};
  M11left=Mmats_v{3};
  M11right=Mmats_v{4};
  %%
  CD=cos(TH-TH0);

  for j=0:N
    jvec0m=(1:Nunc)+j*(2*Nunc);%% a_m row
    jvec1m=(1:Nunc)+(Nunc)+j*(2*Nunc);%% b_m row

    for r=0:N
      jvec0n=(1:Nunc)+r*(2*Nunc);%% a_n column
      jvec1n=(1:Nunc)+(Nunc)+r*(2*Nunc);%% b_n column

      for p=0:N
        %% CALC'S FOR u EQUATIONS:
        %% integral #1:
        %%  < psi_m || U_{n0} || \psi_n >
        KMAT(jvec0m,jvec0n) = KMAT(jvec0m,jvec0n) + ...
          M00left(j+1,p+1)*M00right(p+1,r+1)*dn0U_mn(:,:,p+1);

        %% integral #2:
        %%  - < psi_m || U || \psi_n >
        %%  NB this gives a strong diagonal contribution
        %%     due to the log singularity in the Green's function
        KMAT(jvec0m,jvec1n) = KMAT(jvec0m,jvec1n) + ...
          M01left(j+1,p+1)*M01right(p+1,r+1)*U_mn(:,:,p+1);
      end
    end
  end
end
return


%%
%  [U_nonsing,dn0U,fac_ds0U,U_nonsing_end]=...
%     POOL_hankel_stuff(matrices,vectors,LC,end_points);
%%
%  ip_psi=ip_stuff{1}{3};
%  ip_d1psi=ip_stuff{1}{2};
%  ip_d2psi=ip_stuff{1}{1};
%%
[S0onLC,SonLC]=meshgrid(svec_on_LC,svec_on_LC);
DSonLC=abs(SonLC-S0onLC);
Log0=0*R;
Log0B=0*R;
Log1=0*R;
jnz=find(R);
%%
Log0(jnz)=(1/2/pi)*log( R(jnz)./DSonLC(jnz) );
Log0_mn=ip_d2psi*Log0*ip_d2psi' + log_stuff{1};
%%
Log0B(jnz)=(1/2/pi)*(CD(jnz)-1).*log(DSonLC(jnz));
Log0B_mn=(LC*ip_d1psi)*Log0*(LC*[ip0;ip_d1psi])';
%%
Log1(jnz)=(LC/2/pi)*DSonLC(jnz).*( log(DSonLC(jnz))-1 );
Log1_mn=(LC*ip_d1psi)*Log1*(ip_d2psi)';
%%
Log1plus=(LC/2/pi)*(1-svec_on_LC).*( log(1-svec_on_LC)-1 );
Log1minus=(LC/2/pi)*(1+svec_on_LC).*( log(1+svec_on_LC)-1 );
Log1plus_m=(LC*ip_d1psi)*Log1plus;
Log1minus_m=(LC*ip_d1psi)*Log1minus;

for j=0:N
  jvec0m=(1:Nunc+1)+j*(2*Nunc+1);%% [A,a_m] row
  jj0m=jvec0m(2:end);
  j_Am=jvec0m(1);
  j_a0m=jvec0m(2);
  jvec1m=(1:Nunc)+(Nunc+1)+j*(2*Nunc+1);%% b_m row
%    KMAT_in(jvec0,jvec1) = KMAT_in(jvec0,jvec1) + Mdiag;
  for r=0:N
    jvec0n=(1:Nunc+1)+r*(2*Nunc+1);%% [A,a_n] column
    jj0n=jvec0n(2:end);
    j_An=jvec0n(1);
    j_a0n=jvec0n(2);
    jvec1n=(1:Nunc)+(Nunc+1)+r*(2*Nunc+1);%% b_n column
%      j_An=1+jvec1n(end);%% A col
%      jA(r)=j_An;
    %%
    for p=0:N
      %% CALC'S FOR u EQUATIONS:
      %% integral #1:
      %%  < psi'' || \Uop_-*\pa_{n_0}\Gop_-*\Uop_-^{-1}...
      %%         || \bfu >
      dn0U_mn=(ip_d2psi)*dn0U(:,:,p+1)*(LC*[ip0;ip_d1psi])';
      KMAT(jj0m,jvec0n) = KMAT(jj0m,jvec0n) + ...
                                M00left(j+1,p+1)*M00right(p+1,r+1)*dn0U_mn;

      %% integral #2:
      %%  - < psi'' || \Uop_-*\Gop_-*\Vop_-^{-1}...
      %%           || \bfv >
      %%  NB this gives a strong diagonal contribution
      %%     due to the log singularity in the Green's function
      U_mn=(ip_d2psi)*U_nonsing(:,:,p+1)*(ip_d2psi)' + Log0_mn;%%!!!!!!!!!
      KMAT(jj0m,jvec1n) = KMAT(jj0m,jvec1n) + ...
                                M01left(j+1,p+1)*M01right(p+1,r+1)*U_mn;

      %% CALC'S FOR v EQUATIONS:
      %% integral #1:
      %%  - < LC*psi' || \Vop_-*\pa_n\Gop_-*\Vop_-^{-1}...
      %%             || \bfv >
      dnU_mn=(LC*ip_d1psi)*dn0U(:,:,p+1).'*(ip_d2psi)';
      KMAT(jvec1m,jvec1n) = KMAT(jvec1m,jvec1n) + ...
                                M11left(j+1,p+1)*M11right(p+1,r+1)*dnU_mn;

      %% integral #2:
      %% < LC*psi' || \Vop_-*KK*\Uop_-^{-1}...
      %%            || \bfu >
      %% KK=\cos(\Delta)*\Kop_-^2*( \Gop_- - 1/2/pi*log|(s-s_0)/LC| )
      U2=k_wtr(p+1)^2*CD.*U_nonsing(:,:,p+1);
      U_mn2=(LC*ip_d1psi)*U2*(LC*[ip0;ip_d1psi])';
      KMAT(jvec1m,jvec0n) = KMAT(jvec1m,jvec0n) + ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*U_mn2;

      %% integral #3:
      %% < LC*psi' || \Vop_-*KK*\Uop_-^{-1}...
      %%            || \bfu >
      %% KK=[\cos(\Delta)-1]*\Kop_-^2*( 1/2/pi*log|(s-s_0)/LC| )
      KMAT(jvec1m,jvec0n) = KMAT(jvec1m,jvec0n) + ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  k_wtr(p+1)^2*Log0B_mn;

      %% integral #4:
      %% < LC*psi' || (1/2/pi)*\Vop_-*\Kop_-^2*\log|s-s_0|*\Uop_-^{-1}...
      %%            || \bfu >
      %%  - do this by integrating by parts.
      KMAT(jvec1m,jj0n) = KMAT(jvec1m,jj0n) - ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  k_wtr(p+1)^2*Log1_mn;
      KMAT(jvec1m,j_An) = KMAT(jvec1m,j_An) + ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  k_wtr(p+1)^2*(Log1plus_m-Log1minus_m);
      KMAT(jvec1m,j_a0n) = KMAT(jvec1m,j_a0n) + ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  k_wtr(p+1)^2*Log1plus_m;

      %% integral #5:
      %% < LC*psi' || \pa_s\Vop_-*\Gop_-*\Uop_-^{-1}...
      %%            || \pa_{s_0}\bfu >
      %%  - get rid of \pa_s by integrating by parts:
      %%
      %%  NB this gives a strong diagonal contribution
      %%     due to the log singularity in the Green's function
      KMAT(jvec1m,jj0n) = KMAT(jvec1m,jj0n) - ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  U_mn;
      %%  NB integration by parts results in an extra term for m=0:
      U_m_end=log_stuff{2}+U_nonsing_end.'*ip_d2psi';
      KMAT(jvec1m(1),jj0n) = KMAT(jvec1m(1),jj0n) + ...
                                M10left(j+1,p+1)*M10right(p+1,r+1)*...
                                  U_m_end;
    end
  end
end



if ~USE_XTRA%% don't need the xtra constant A
  Kmat_u(:,jA)=[];
end