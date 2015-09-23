function run_SlepianD(aus,testEnergy)
%run_SlepianD: compute and plot Slepian eigenfunctions on the sphere (8.27) and (8.29)

if nargin<1
	aus=1; % default is to do Australia
end
if nargin<2
	testEnergy=0;
end

base=userpath;
base(end)='/'; % ~/Documents/MATLAB
frames_folder=[base 'frames/'];

intginc=0.1; % fine grid for integration (degrees)
medinc=0.5; % medium grid for smooth plotting (degrees)
plotinc=2.0; % coarse grid for grid plotting (degrees)

fprintf('\n@@ Fine Grid: %.2f (degrees)\n',intginc)
fprintf('@@ Medium Grid: %.2f (degrees)\n',medinc)
fprintf('@@ Coarse Grid: %.2f (degrees)\n',plotinc)

close all

if aus % region is Australia including Tasmania
	[tv_aus,pv_aus,mask_aus,tR,pR]=ausRegion(medinc,1);
	basename='aus';
	% tv_aus,pv_aus are the bounding vectors of the rectangle containing
	% Australia. mask_aus is the 0-1 matrix indicating Australia within the
	% rectangle. tR,pR define the coastline with NaN separating subregions.
else % any number of non-intersecting regular subregions
	% subregion 1
	tr(1,:)=[30 60];
	pr(1,:)=[-70 -10];
	% subregion 2
	tr(2,:)=[80 110];
	pr(2,:)=[-70 -10];
	% subregion 3
	tr(3,:)=[20 70];
	pr(3,:)=[0 40];
	basename='reg';
end

for L_max=3:3 % range of L_max
	fprintf('\n@@ L_max: %d\n',L_max)
	% populate the N_tot x N_tot Hermitian D matrix (8.27)
	N_tot=(L_max+1)^2;
	if aus % irregular non-simply-connected (uses mask)
		D=SlepianDH(L_max,tv_aus,pv_aus,mask_aus);
	else
		D=zeros(N_tot,N_tot); % allocate and initialize to zero
		for r=1:size(tr,1) % loop over subregions
			tv_intg=(tr(r,1):intginc:tr(r,2))*pi/180;
			pv_intg=(pr(r,1):intginc:pr(r,2))*pi/180;
			D=D+SlepianDH(L_max,tv_intg,pv_intg); % accumulate
			fprintf('@@ Completed D for Region: %d\n',r)
		end
	end
	fprintf('@@ Size of D matrix: %dx%d\n', size(D))

	% eigindex'th eigenvector w with eigenvalue lambda in spectral domain
	[V,lamD]=eig(D); % eigen-structure
	% stem(flip(diag(lamD))) % plot the eigenvalues

	for eigindex=0:N_tot-1 %0:3%N_tot-1 % 0:N_tot-1 eigenvalue index 0,1,2 in descending energy order
		fprintf('\n@@ Eigenindex: %d\n',eigindex)
		lambda=min(max(lamD(end-eigindex,end-eigindex),0.0),1.0);
		fprintf('@@ Eigenvalue %d: %8.6f\n',eigindex,lambda)

		%% spectral Slepian eigenvector
		w=V(:,end-eigindex); % eigenvector - SHT of Slepian function
		w=w/norm(w); % normalize eigenvector (actually this is superfluous)

		if testEnergy~=0 % test spatial energy of slepian function is unity
			spectralEnergy=norm(w)^2; % will be 1 but just check
			fprintf('@@ Spectral Energy: %8.6f\n',spectralEnergy)

			% eigindex'th Slepian eigenfunction in spatial domain (0'th is dominant)
			tv_intg=(0:intginc:180)*pi/180;
			pv_intg=(0:intginc:360)*pi/180;

			%% spatial Slepian eigenfunction [slepian,theta,phi]
			[slepian,~,~]=spatial(w,tv_intg,pv_intg); % ISHT; theta,phi unused

			spatialEnergy=trapSphereMaskedR(abs(slepian).^2,tv_intg,pv_intg);
			fprintf('@@ Spatial Energy: %8.6f\n',spatialEnergy)
		end

		%% plot eigenfunction in spatial domain
		close % prepare fresh figure
		bump_height=0.15;
		ref_sphere=1.0;
		plottype=1;
		tv_plot=(0:medinc:180)*pi/180;
		pv_plot=(0:medinc:360)*pi/180;
		[slepian_plot,theta_plot,phi_plot]=spatial(w,tv_plot,pv_plot,1);
		maxSlepian=max(abs(slepian_plot(:))); % for plot normalization
		fprintf('@@ Max Slepian: %8.6f\n',maxSlepian)
		s=spatial_plot(slepian_plot/maxSlepian,theta_plot,phi_plot, ...
			bump_height,ref_sphere,plottype);
		s.EdgeColor='none'; % no lines
		s.FaceAlpha=0.8;
		colormap('cool') % or copper
		hold on

		if aus==1
			% fill Australia area on zero/reference surface
			aus_coast=ref_sphere*[sin(tR).*cos(pR); sin(tR).*sin(pR); cos(tR)]/ ...
				(ref_sphere+bump_height);
			idx = any([~isnan(tR);~isnan(tR)],1);
			fill3(aus_coast(1,idx),aus_coast(2,idx),aus_coast(3,idx), ...
				'yellow','EdgeColor','None');
			% draw Australia coastline on Slepian surface
			rR=interpn(theta_plot,phi_plot,slepian_plot,tR,pR); % very slick
			F=spatialMap(rR/maxSlepian,plottype);
			rad=abs(ref_sphere + bump_height*F)/(ref_sphere+bump_height);
			aus_coast=[rad.*sin(tR).*cos(pR); rad.*sin(tR).*sin(pR); rad.*cos(tR)];
			plot3(aus_coast(1,:),aus_coast(2,:),aus_coast(3,:),'yellow','LineWidth',2);
			view(-100,-30) % set viewpoint
			set(gca,'CameraViewAngle',9) % zoom into scene
		else
			% Compute energy in sub-regions
			lambda_est=0;
			for r=1:size(tr,1) % loop over subregions (fine grid)
				tv_intg=(tr(r,1):intginc:tr(r,2))*pi/180;
				pv_intg=(pr(r,1):intginc:pr(r,2))*pi/180;
				reg=spatial(w,tv_intg,pv_intg,0);
				lambda_est=lambda_est+trapSphereMaskedR(abs(reg).^2,tv_intg,pv_intg);
			end
			fprintf('@@ Region Energy: %8.6f (eigenvalue %8.6f)\n',lambda_est,lambda)

			% Plot sub-regions
			for r=1:size(tr,1) % loop over subregions (coarse grid)
				ttcr=(tr(r,1):plotinc:tr(r,2))*pi/180;
				ppcr=(pr(r,1):plotinc:pr(r,2))*pi/180;
				[r_reg,theta_reg,phi_reg]=spatial(w,ttcr,ppcr,0);
				s=spatial_plot(1.01*r_reg/maxSlepian,theta_reg,phi_reg);
				s.EdgeColor='yellow'; % no lines
			end
		end

		% Annotate plot
		llabel=sprintf('$L_{\\mathrm{max}}=%d$\n$\\lambda_{%d}=%8.6f$', ...
			L_max,eigindex,lambda);
		delete(findall(gcf,'Tag','myLabel'));
		a=annotation('textbox',[0.695,0.895,0.28,0.1],'String',llabel);
		a.Interpreter='latex';
		a.FontSize=18;
		a.LineStyle='none';
		a.Tag='myLabel';
		hold off

		% output to png file (directory needs to exist)
		outname=sprintf('%s_%04d_%04d',[frames_folder basename],L_max,eigindex);
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
		saveas(gcf,outname,'png')
   end

	%% render eigenfunction movie on osx; needs avconv - get via brew install libav
	if ~system('which avconv >/dev/null')
		renderMovWithAvconv=sprintf(...
			'avconv -framerate 2 -y -v quiet -f image2 -i %s_%04d_%%04d.png %sm-%04d.mov',...
			[frames_folder basename],L_max,[frames_folder basename],L_max);
		system(renderMovWithAvconv); % make compressed mov
      cleanup=sprintf('rm %s_%04d*.png',[frames_folder basename],L_max);
		system(cleanup); % delete frames
	end
end
