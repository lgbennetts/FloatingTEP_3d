function mx = Max_Floe_new(edge_vars, res, conts, Psi, Psi_dn,...
    kk, mat_C, Vert_Dim, Nmodes, Nint, fig)

%% - 04.03.10: rewritten by LGB -> extract sings in Green's fns - %%
%% - 01.12.08: written by LGB - %%
% Nint = 1*300;

NN = -Nmodes:Nmodes;

xmat = zeros(conts(1)+1+2, res+1);
ymat = zeros(conts(1)+1+2, res+1);
% wmat = zeros(conts(1)+1+2, res+1);
whatmat = zeros(conts(1)+1+2, res+1);

cont_vec = linspace(0, 1, conts(1)+1+2);
% cont_vec = sqrt(cont_vec);

%% - mod of Tim's stuff for generating curve - use for ints
 [tt,ww]=OP_numint_legendre(Nint);%% use Gauss-legendre points
                          %% and weights as everything smooth.
[xyvecs_bdy,ds_dt_bdy,th_vec_bdy,~,~,~,~,LC_bdy,svec_bdy]=CRK_get_rsdtheta(edge_vars,tt);
xvec=xyvecs_bdy(1,:).';
yvec=xyvecs_bdy(2,:).';
[~,X0_bdy]=meshgrid(xvec,xvec);
[~,Y0_bdy]=meshgrid(yvec,yvec);
[~,TH0_bdy]=meshgrid(th_vec_bdy,th_vec_bdy); 
[~,S0_bdy]=meshgrid(svec_bdy,svec_bdy);

clear xvec yvec X Y

ww_bdy=ww.*ds_dt_bdy/LC_bdy;
[IP_bdy,~,~]=GEN_inprod_exp(svec_bdy/LC_bdy,ww_bdy,Nmodes);

lam_n = -1i*pi*NN/LC_bdy; % - size 2N+1-sq (for diffing the exp basis fns)

IPds0_bdy = diag(lam_n)*IP_bdy; 

%% - Inner Floe - %
%% -- now for the fn evals & fill 1st entries with boundary vals -- %
 
res_vec = linspace(-1,1,res+1);
[dum_xyvecs,~,~,~,~,~,~,dum_LC,dum_svec]=CRK_get_rsdtheta(edge_vars,res_vec.');
xmat(conts(1)+1+2,:)=dum_xyvecs(1,:).';
ymat(conts(1)+1+2,:)=dum_xyvecs(2,:).';

dum_Psi =  (1/2/dum_LC)*exp(1i*dum_svec*NN*pi/dum_LC)*Psi.'; dum_Psi = dum_Psi.'; % - soln on the bdy (already known)
% wmat(1, :) = dum_Psi(Vert_Dim+1,:);  clear dum_Psi
whatmat(1, :) = dum_Psi(Vert_Dim+2,:);  clear dum_Psi

%% -- loop points on inner contours

edge_vars2 = edge_vars; % -> for calc'ing inner contours

invCPsi = mat_C\Psi; invCPsi_dn = mat_C\Psi_dn;

for  loop_c = conts(1)+1+1:-1:2 % -> run from outside inwards
    
    invCPsi_in = 0*Psi;
    
    clear dum_Psi
    edge_vars2{3}{1} = cont_vec(loop_c)*edge_vars{3}{1};
    
    [xyvecs_in,ds_dt_in,~,~,~,~,~,LC_in,svec_in]=CRK_get_rsdtheta(edge_vars2,tt);

    xvec=xyvecs_in(1,:).'; yvec=xyvecs_in(2,:).';
    [X_in,~]=meshgrid(xvec,xvec); [Y_in,~]=meshgrid(yvec,yvec);
    clear xvec yvec
    [S_in,~]=meshgrid(svec_in,svec_in);

    [R,varTH]=Polars_cts(X0_bdy-X_in,Y0_bdy-Y_in);
    
    ww_in=ww.*ds_dt_in/LC_in;
    [IP_in,~,~]=GEN_inprod_exp(svec_in/LC_in,ww_in,Nmodes);
    
    % - The singular bit of H0: (2i/pi)log(R) R^2=(x-x0)^2 + (y-y0)^2
    
    LogOn = 1; % - INNER!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    if LogOn == 1
        % - v1 - %
        %LogTerm = log(abs(1-(LC_in/LC_bdy)*exp(1i*pi*(S0_bdy/LC_bdy - S_in/LC_in)))) - 0.5*log(pi^2/LC_in/LC_bdy);
        % ------ %
        % - v2 - %
        LogTerm = log(abs(1-(LC_in/LC_bdy)*exp(1i*pi*(S0_bdy/LC_bdy - S_in/LC_in)))) - log(pi/LC_bdy) - log(LC_in/LC_bdy);
        % ------ %
        
        LogTerm = (2i/pi)*LogTerm;
        
        SingTerm = (2i/pi)*sin(TH0_bdy-varTH)./R;

        dum_f = -1i*exp(1i*TH0_bdy).*(1-(LC_in/LC_bdy)*exp(1i*pi*(S0_bdy/LC_bdy - S_in/LC_in))); % - nb n0=-1i*exp(1i*TH0_bdy)
        [~,THf]=Polars_cts(real(dum_f),imag(dum_f));
        resTH = (2i/pi)*(varTH+THf);
        
        % ---------------------------- %
        
        dum_vec = (1i/pi)*((LC_in/LC_bdy).^abs(NN))./abs(NN);
        % - v1 - %
        %dum_vec(Nmodes+1) = 1i*log(pi^2/LC_in/LC_bdy)/pi;
        % ------ %
        % - v2 - %
        dum_vec(Nmodes+1) = 2i*(log(pi/LC_bdy)+log(LC_in/LC_bdy))/pi;
        % ------ %
        dum_vec = dum_vec/4i;
            
        IBP_mat = zeros(length(NN));
        IBP_mat(:,Nmodes+1) = (2i/LC_bdy)*(-1).^NN; IBP_mat = IBP_mat/4i;
        
        int_resTH = conj(IPds0_bdy)*resTH*IP_in.'; int_resTH = int_resTH/4i;
        
        int_TH0 = conj(IPds0_bdy)*unwrap(unwrap(TH0_bdy,1),2)*IP_in.'; int_TH0 = int_TH0/2/pi;
        
        dum_vecdn = -lam_n.*((LC_in/LC_bdy).^abs(NN))./NN;
        dum_vecdn(Nmodes+1) = 0; dum_vecdn = dum_vecdn/4i/pi;
            
    end
    
    for loop_Dim=1:Vert_Dim+2
        
        if LogOn == 1
            H0 = besselh(0, kk(loop_Dim)*R)-LogTerm;
            H0_dn = -kk(loop_Dim)*sin(TH0_bdy-varTH).*besselh(1, kk(loop_Dim)*R) ...
                - SingTerm; 
        elseif LogOn==0
            H0 = besselh(0, kk(loop_Dim)*R);
            H0_dn = -kk(loop_Dim)*sin(TH0_bdy-varTH).*besselh(1, kk(loop_Dim)*R); 
        end
        
        H0 = H0/4i; H0_dn = H0_dn/4i;
        
        H0I = conj(IP_bdy)*H0*IP_in.'; H0nI = conj(IP_bdy)*H0_dn*IP_in.';
        
        clear H0 H0_dn
        
        % - correct for log
        
        if LogOn == 1
            
            H0I = H0I - diag(dum_vec); 
            
            % ----------------- %
            
            H0nI = H0nI + IBP_mat - int_resTH + int_TH0 + diag(dum_vecdn);
            
        end
        
        invCPsi_in(loop_Dim,:) = (transpose(H0nI)*invCPsi(loop_Dim,:).' - transpose(H0I)*invCPsi_dn(loop_Dim,:).').';
        
    end
    
    clear LogTerm SingTerm dum_vec dum_vecdn int_TH0 int_LogTH IBP_mat
            
    [dum_xyvecs,~,~,~,~,~,~,dum_LC,dum_svec]=CRK_get_rsdtheta(edge_vars2,res_vec.');
    xmat(loop_c,:)=dum_xyvecs(1,:).';
    ymat(loop_c,:)=dum_xyvecs(2,:).';

    dum_Psi(1:Vert_Dim+2,:) = (mat_C*invCPsi_in)*exp(1i*dum_svec*NN*pi/dum_LC).';

    %wmat(loop_c, :) = dum_Psi(Vert_Dim+1,:);
    whatmat(loop_c, :) = dum_Psi(Vert_Dim+1,:);
    
    clear H0I H0nI invCPsi_in xyvecs_in ds_dt_in LC_in svec_in X_in Y_in S_in R varTH
    
end

clear dum_Psi

edge_vars2{3}{1} = cont_vec(1)*edge_vars{3}{1};
[dum_xyvecs,~,~,~,~,~,~,~,~]=CRK_get_rsdtheta(edge_vars2,res_vec.');
xmat(1,:)=dum_xyvecs(1,:).'; ymat(1,:)=dum_xyvecs(2,:).';

xx = xyvecs_bdy(1,:).'-xmat(1,1); yy = xyvecs_bdy(2,:).'-ymat(1,1); zz = xx + 1i*yy; % -> origin (treat once only)

H0 = besselh(0, abs(zz)*kk);
H0_dn = -((1./abs(zz)).*( sin(th_vec_bdy).*xx - cos(th_vec_bdy).*yy )*kk).*...
      ...
      besselh(1, abs(zz)*kk);
  
H0 = H0/4i; H0_dn = H0_dn/4i;
  
H0I = conj(IP_bdy)*H0; H0nI = conj(IP_bdy)*H0_dn;
  
dum_Psi =  mat_C*(transpose(H0nI).*(invCPsi) ...
      - transpose(H0I).*(invCPsi_dn))*ones(2*Nmodes+1,1); 
% wmat(1, :) = dum_Psi(Vert_Dim+1);
whatmat(1, :) = dum_Psi(Vert_Dim+2);

clear H0 H0_dn HOI H0nI xx yy zz

if 0
    %% - plot - %%
    figure(fig); hold on
    
    %% -- Pool -- %%
    
    % [c, h] = contour(xmat, ymat, real(dispmat),'k'); clabel(c,h);
    
    surface(xmat, ymat, real(wmat))
    shading interp
    colormap gray
    % set(gcf,'colormap', get(gcf,'colormap')/2 + 0.25);
    % contour3(xmat, ymat, real(dispmat),'k')
    
    %% - Ice -- %%
    
    % [c, h] = contour(xmatI, ymatI, real(wmat),'k','linestyle',':'); clabel(c,h);
    
    surface(xmatI, ymatI, real(dispmat))
    shading interp
    
    set(gca, 'xlim', 40*[-1, 1]);
    set(gca, 'ylim', 40*[-1, 1]);
    
    xlabel('x'); ylabel('y');
    
end

mx = max(max(abs(real(whatmat))));