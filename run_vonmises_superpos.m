function run_vonmises_superpos
%run_Rotsymtest render a superposition of symmetric functions

close all
clear ll w_eta w_mu w

L = 30; % maximum (inclusive) degree

%% mu - arbitrary angle component

theta = [0,0.55]*pi;
phi = [0,1.5]*pi;
kappa = [30,90];
alpha = [1,1];


%% SH coefficients of superposition
w = 0;
for i=1:length(alpha)
    w = w + alpha(i) * vonmises(theta(i),phi(i),kappa(i),L); %superposition
end

%% determine spatial samples on mesh
deginc=3.0;
tv=(0:deginc:180)*pi/180;
pv=(-30:deginc:330)*pi/180;
[w_spatial,theta,phi]=ishtRectGrid(w,tv,pv,true);

%% do plot
doThePlot(w_spatial,theta,phi,'rotationally symmetric function plot');


function doThePlot(Ylm,theta,phi,name)
%doThePlot - convenience function to plot the spatial form
maxY=max(abs(Ylm(:)));
Ylm=Ylm/maxY; % normalize
figure
bump_height=0.8; ref_sphere=1.0;
plottype=2; % real
s=spatialPlot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='none'; % 'none'
fig=gcf;
fig.Name=name;
end

function [w_mu] = vonmises(theta_mu,phi_mu,kappa_mu,L)
% range of degree
N_tot=(L+1)^2; % max number of SH coefficients
ll=0:L; % index 1,2,...,L+1 convenience vector    
w_mu=zeros(N_tot,1); % 1D coeffs initially filled with zeros
lambda_mu=exp(-ll.*(ll+1)/(2*kappa_mu)); % Gauss-Weierstrass kernel
% use formula to determine the SH coefficients of mu componentdf
for l=0:L
	for m=-l:l
		n=l*(l+1)+m; % (7.39)
		w_mu(n+1)=lambda_mu(l+1)*conj(sphHarm(l,m,theta_mu,phi_mu)); % scalar angles here
	end
end
end
end
