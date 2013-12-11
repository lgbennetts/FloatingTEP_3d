function [rvecs,Svec,th_vec, dtheta_ds,dS_dt,d2theta_ds2,...
   d2r_ds2,LC,tvec]=CRK_get_rsdtheta_NR(crk_vars,svec_on_LC);

fxn_handle=crk_vars{1};%%which routine to use for crack shape
crk_prams=crk_vars{2};%%parameters needed for given shape
srt=crk_vars{3};%%scaling,rotation,translations

%% calc length of curve:
s_target=0;
LC=feval(@arc_length,1,fxn_handle,crk_prams,srt,s_target)/2;

%% now find which values of t give desired value of s,
%% and thus find the required quadrature points.
tvec=0*svec_on_LC;
interval=[-1 1];
guess=svec_on_LC(end);

for j=length(svec_on_LC):-1:1
  s_target=(1+svec_on_LC(j))*LC;
  [tvec(j),flag]=GEN_findroot_NRsafe(...
     @arc_length,interval,guess,fxn_handle,crk_prams,srt,s_target);
  %flag
  interval=[-1 tvec(j)];
  guess=-1+.99*(1+tvec(j));
end



%for j=1:length(tvec)
%  tst(j,:)=[svec_on_LC(j)*LC,...
%    feval(@arc_length,tvec(j),fxn_handle,crk_prams,srt,0)-LC];
%end
%tst,LC,return

%% calc r,dr,d2r for basic function
[dr,rvecs,d2r,d3r]=feval(fxn_handle,tvec,crk_prams);

%% then allow for scaling, rotations and translations
scaler=srt{1};
M_sc=diag(scaler);
th=pi*srt{2}/180;
M_rot=[cos(th),-sin(th);sin(th),cos(th)];
v_trans=srt{3};
rvecs=diag(v_trans)*ones(size(rvecs))+M_rot*M_sc*rvecs;
%%
dr=M_rot*M_sc*dr;
dx=dr(1,:)';
dy=dr(2,:)';
[Svec,th_vec]=GEN_polar_coords(dx,dy);
%[Svec,th_vec]
%tangvecs=[ dx./Svec dy./Svec ]';
%%
d2r=M_rot*M_sc*d2r;
d2x=d2r(1,:)';
d2y=d2r(2,:)';
dtheta_ds=( d2y.*dx - d2x.*dy )./Svec.^3;
dS_dt=(d2x.*dx+d2y.*dy)./Svec;
%%
d3r=M_sc*d3r;
d3x=d3r(1,:).';
d3y=d3r(2,:).';
d2theta_dt2=(d3y.*dx-d3x.*dy)./Svec.^2-2*dS_dt.*dtheta_ds;
d2theta_ds2=(d2theta_dt2-dS_dt.*dtheta_ds)./Svec.^2;
%%
d2r_ds2=[(d2x.*Svec-dS_dt.*dx)./Svec.^3,...
                (d2y.*Svec-dS_dt.*dy)./Svec.^3];

%% integrate wrt tau=s/LC
%%  => ds/d(tau)=LC & d^2s/dtau^2=d(LC)/d(tau)=0:
dS_dt=0*svec_on_LC;
Svec=LC+dS_dt;

function [f,f_on_df]=...
           arc_length(t,fxn_handle,crk_prams,srt,s_target)
%% f=s-s_target

scaler=srt{1};%% don't need to worry about rotations
              %% & translations when calc'ing arc length:
M_sc=diag(scaler);
dr=M_sc*feval(fxn_handle,t,crk_prams);
ds=sqrt( dr(1,:).^2+dr(2,:).^2 );
%%
Ngl=0;
tol=1e-12;
err=2*tol;
s0=0;
while err>tol
  Ngl=Ngl+50;
  [tgl,wgl]=OP_numint_legendre(Ngl,[-1 t]);
  dr=M_sc*feval(fxn_handle,tgl,crk_prams);
  ds_=sqrt( dr(1,:).^2+dr(2,:).^2 );
  s=ds_*wgl;
  %%
  f=s-s_target;
  f_on_df=f/ds;
  err=abs(1-s0/s);
  s0=s;
end