
% load timecourse file
[file,folder] = uigetfile('*_timecourse.mat');
load([folder file])
% also load the descriptor
load([folder strrep(file,'_timecourse','_descriptor')])

n_ex = numel(descriptor.excitationWavelength);
n_em = numel(descriptor.wl_emission);
n_cl = numel(descriptor.Chloride);
n_ph = numel(descriptor.pH);
assert(n_cl==n_ph,'cl values are more/less than ph values!')
n_cond = n_cl;
n_roi = numel(labels);
vals = reshape(timecourse,n_ex,n_em,n_cond,n_roi);

% background
last_roi_bkg = true;
if last_roi_bkg
    vals_noBkg = vals-vals(:,:,:,end);
    vals_noBkg = vals_noBkg(:,:,:,1:end-1);
    n_roi = n_roi-1;
else
    vals_noBkg = vals;
end

wl_em = descriptor.wl_emission;
wl_ex = descriptor.excitationWavelength;

% compute average on rois
avg = squeeze(mean(vals(:,:,:,:),4)); % avg is excitation*emission*condition
normalize = true;
% normalize = false;
if normalize
    avg = avg./mean( avg(:,wl_em>=650,:) , 2);
end
ymax = max(avg,[],'all');

%% dividi per clori
figure
cols = turbo(n_ex);
sbpn = numSubplots(n_cond+1);
for i=1:n_cond
    subplot(sbpn(1),sbpn(2),i)
    hold on
    avg_now = avg(:,:,i); % avg_now is excitation*emission
    for j=1:n_ex
        plot(wl_em, avg_now(j,:), 'color',cols(j,:),'linewidth',1.5)
    end
    title(sprintf('pH = %.2f, [Cl] = %.1f',descriptor.pH(i),descriptor.Chloride(i)))
    ylim([0 ymax])
    xlabel('Emission wl (nm)')
    ylabel('Intensity')
end
subplot(sbpn(1),sbpn(2),n_cond+1)
colormap(cols)
cb=colorbar('east');
set(gca,'CLim',[0.5 n_ex+0.5]);
cb.Ticks = 1:n_ex;
cb.TickLabels = num2str(descriptor.excitationWavelength');
cb.Label.String = 'Ex wl (nm)';
cb.Label.FontSize = 16;

%% dividi per eccitaz
figure
cols = turbo(n_cond);
sbpn = numSubplots(n_ex+1);
for i=1:n_ex
    subplot(sbpn(1),sbpn(2),i)
    hold on
    avg_now = squeeze(avg(i,:,:)); % avg_now is emission*condition
    for j=1:n_cond
        plot(wl_em, avg_now(:,j), 'color',cols(j,:),'linewidth',1.5)
    end
    title(sprintf('Exc = %.1f nm', descriptor.excitationWavelength(i)))
    ylim([0 ymax])
    xlabel('Emission wl (nm)')
    ylabel('Intensity')
end
subplot(sbpn(1),sbpn(2),n_ex+1)
colormap(cols)
cb=colorbar('east');
set(gca,'CLim',[0.5 n_cond+0.5]);
cb.Ticks = 1:n_cond;
cb.TickLabels = num2str(descriptor.Chloride');
cb.Label.String = '[Cl] (Mm)';
cb.Label.FontSize = 16;
set(gca,'visible','off')

%% superf

% visualize = 'image';
visualize = 'surf';
sbpn = numSubplots(n_cond);
switch visualize
    case 'surf'
        [EX,EM] = meshgrid(wl_ex, wl_em);
        figure
        for ic=1:n_cond
            subplot(sbpn(1),sbpn(2),ic)
            surf(EX,EM,squeeze(avg(:,:,ic))','EdgeColor','k','FaceColor','interp')
            xlabel('Ex wl (nm)')
            ylabel('Em wl (nm)')
            zlabel('Norm. intensity')
            colorbar
            if normalize
                set(gca,'clim',[0 6],'View',[0,90])
            else
                set(gca,'clim',[0 26000],'View',[0,90])
            end
            title(sprintf('pH = %.2f, [Cl] = %.1f',descriptor.pH(ic),descriptor.Chloride(ic)))
        end
        
        
    case 'image'
        figure
        subplot(1,2,1)
        imagesc(wl_ex,wl_em,squeeze(avg(:,:,1))')
        xlabel('Ex wl (nm)')
        ylabel('Em wl (nm)')
        colorbar
        if normalize
            set(gca,'clim',[0 6],'ydir','normal')
        else
            set(gca,'clim',[0 26000],'ydir','normal')
        end
        title('Cloro 0')
        
        subplot(1,2,2)
        imagesc(wl_ex,wl_em,squeeze(avg(:,:,end))')
        xlabel('Ex wl (nm)')
        ylabel('Em wl (nm)')
        colorbar
        if normalize
            set(gca,'clim',[0 6],'ydir','normal')
        else
            set(gca,'clim',[0 26000],'ydir','normal')
        end
        title('Cloro 40')
end

%% spettri emissione e kd

Rem = avg(1,:,end);
Gem = max(avg(1,:,2)-Rem,0);

% normalize to peak
Rem = Rem./max(Rem);
Gem = Gem./max(Gem);

figure
plot(wl_em,Gem,'g')
hold on
plot(wl_em,Rem,'r')

coef0 = [0.5 0.5];
coef = nan(numel(wl_ex),2,n_cond);

for i_ex = 1:numel(wl_ex)
    Iem = squeeze(avg(i_ex,:,:));
    for i=1:n_cond
        resid = @(x) ( mean( (Iem(:,i) - x(1).*Gem' - x(2).*Rem').^2, 'all', 'omitnan'));
        coef(i_ex,:,i) = fminsearch(resid,coef0);
    end
    figure('Name',sprintf('Cl fit, exc = %i nm',wl_ex(i)))
    for i=1:n_cond
        subplot(1,n_cond,i)
        plot(wl_em,Iem(:,i),'ok')
        hold on
        plot(wl_em,coef(i_ex,1,i).*Gem,'g')
        plot(wl_em,coef(i_ex,2,i).*Rem,'r')
        plot(wl_em,coef(i_ex,1,i).*Gem+coef(i_ex,2,i).*Rem,'m')
        title(sprintf('[Cl] = %.1f',descriptor.Chloride(i)))
    end
end

cl = descriptor.Chloride;

figure
z = squeeze( coef(1,1,:) ./ ( coef(1,1,:)+coef(1,2,:) ) )';
plot(cl, z, 'o-')
xlabel('[Cl] (mM)')

rcl = @(cl,r0,kd) (r0 .* 1 ./ (1 + cl./kd));
resid2 = @(x) ( mean( (z - rcl(cl,x(1),x(2)) ).^2 ) );
par0 = [0.5 5];
par = fminsearch(resid2,par0);
hold on
plot(cl(1):1:cl(end), rcl(cl(1):1:cl(end),par(1),par(2)), 'r')


%% excitation spectra

first_wl = 520;
em2 = wl_em( wl_em>=first_wl); 
% BP_G = [492 562];
% BP_R = [572 642];
% alpha = 0.1637; % red into green
% beta = 0.1108; % green into red
BP_G = [520 540];
BP_R = [600 680];
alpha = 0;
beta = 0;
G = squeeze( sum( avg(:,wl_em>=BP_G(1) & wl_em>=first_wl & wl_em<=BP_G(2),:), 2) );
R = squeeze( sum( avg(:,wl_em>=BP_R(1) & wl_em>=first_wl & wl_em<=BP_R(2),:), 2) );

Gstar = G-alpha.*R; % exc*chloride
Rstar = R-beta.*G; % exc*chloride

cols = turbo(n_cond);
wl_ex = descriptor.excitationWavelength;

figure
subplot(3,1,1)
hold on
for i=1:n_cond
    plot(wl_ex,Gstar(:,i),'color',cols(i,:),'linewidth',2)
end
title('Green')

subplot(3,1,2)
hold on
for i=1:n_cond
    plot(wl_ex,Rstar(:,i),'color',cols(i,:),'linewidth',2)
end
title('Red')

subplot(3,1,3)
hold on
for i=1:n_cond
    plot(wl_ex,Gstar(:,i)./Rstar(:,i),'color',cols(i,:),'linewidth',2)
end
title('Green/Red')
