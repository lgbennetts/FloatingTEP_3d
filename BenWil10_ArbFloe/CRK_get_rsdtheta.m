function [rvecs,Svec,th_vec, dtheta_ds,dS_dt,d2theta_ds2,...
   d2r_ds2,LC,svec]=CRK_get_rsdtheta(crk_vars,tt);

fxn_handle=crk_vars{1};%%which routine to use for crack shape
crk_prams=crk_vars{2};%%parameters needed for given shape
srt=crk_vars{3};%%scaling,rotation,translations

%% calc r,dr,d2r for basic function
[dr,rvecs,d2r,d3r]=feval(fxn_handle,tt,crk_prams);

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
% % - LB (26.1.09) - %
% [dum,th_vec]=GEN_polar_coords(rvecs(1,:).',rvecs(2,:).');
%[Svec,th_vec]
%tangvecs=[ dx./Svec dy./Svec ]';
%%
d2r=M_rot*M_sc*d2r;
d2x=d2r(1,:)';
d2y=d2r(2,:)';
dtheta_ds=( d2y.*dx - d2x.*dy )./Svec.^3;
% - LB (27.08.09) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! - %
% dtheta_ds=( d2y.*dx - d2x.*dy )./Svec.^2; 
% dtheta_ds=dtheta_ds/LC;
% ----------------- %
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
%%
if nargout>=8
  %% GET THE CRACK LENGTH:
  Nint=0; Nint_tol=2000;
  LC=0;
  tol=1e-12;
  dL=1;
  while and(dL>tol,Nint<=Nint_tol)
    Nint=Nint+50;
    LC0=LC;
    %%
    [uu,ww]=OP_numint_legendre(Nint);
    dr=M_rot*M_sc*feval(fxn_handle,uu,crk_prams);
    dx=dr(1,:)';
    dy=dr(2,:)';
    S=abs(dx+1i*dy);
    %%
    LC=sum(ww.*S)/2;
    dL=abs(1-LC0/LC);
  end
  if Nint>Nint_tol
   disp(['warning: boundary length not converged: ' num2str(dL)])
  end
end

if nargout==9
  %% NOW EVALUATE THE ARC-LENGTH AT THE QUAD POINTS:
  svec=0*tt;
  for j=1:length(tt)
    t=tt(j);
    tau=-1+(uu+1)*(t+1)/2;
    dr=M_rot*M_sc*feval(fxn_handle,tau,crk_prams);
    dx=dr(1,:)';
    dy=dr(2,:)';
    S=abs(dx+i*dy);
    %%
    svec(j)=(t+1)/2*sum(ww.*S) - LC;
  end
end
