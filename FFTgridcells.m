clar
%addpath(genpath('C:\Users\seth.koenig\Documents\MATLAB\generic'))

%load(['R:\Buffalo Lab\Nathan\backup\RateMaps\' respcell{k,1} 'ratemapFinal2.mat'])
%load('R:\Buffalo Lab\Nathan\backup\RateMaps\MP111213 sig007cratemapFinal2.mat')
%load('R:\Buffalo Lab\Nathan\backup\RateMaps\MP111213 sig004aratemapFinal2.mat')

load('R:\Buffalo Lab\Nathan\backup\RateMaps\TT091218 sig003aratemapFinal2.mat')

numreps = 10;
%%
% Based on rmap: compute ratemap and autocorr input spike train, position values, and options by N Killian 110905

o = setdefaults(o,'getTAC',0,'useNaNfilt',0,'imfiltType',0,'useoldtimesforth',0,'mosersmth',0,'nbins',50,...
    'binningtype','specrangeExact','sqr',ceil(max(abs(Yd))),'degperbin',0.5,'dogauss',1,'gsize',[13 13],'gsigma',3,...
    'spkth',0,'sth',1,'smtht',1,'timth',1,'tth',15,'dsfs',1e3,'mincorrpixels',15,'usexcorr2',0,...
    'keepSinglePixels',0,'Nrbin',0);

rnanind = [];
H = fspecial('gaussian',o.gsize,o.gsigma);

S = Sd;
X = Xd;
Y = Yd;

% do binning then smooth spikes and time, this is the Killian method
% it's computationally much faster than the Moser-style methods,
% but smooths over the kernel range after binning so is potentially less precise
disp('using NK method...')

[h t r xb yb lb rb x0 y0 z0 origind] = hist3gc(X,Y,S,o.nbins,o.nbins,o.binningtype,o.sqr,o.degperbin);

%used for saving original data
h0 = h;
t0 = t;
r0 = r;

% Smooth Time and The Spikes Separately
h2 = imfilter(h,H,o.imfiltType);
t2 = imfilter(t,H,o.imfiltType);

shortind = t2<o.tth; % USE SMOOTHED TIMES FOR THRESHOLDING
h2(shortind) = nan;
t2(shortind) = nan;

r = h2./t2;
I2 = r;

numrow = size(t0,1);
numcol = size(t0,2);

shuffled_rs = cell(1,numreps);
for reps = 1:numreps
    
    shuffind = randperm(numrow*numcol);
    
    sh0 = reshape(h0(shuffind),numrow,numcol);
    st0 = reshape(t0(shuffind),numrow,numcol);
    
    sh2 = imfilter(sh0,H,o.imfiltType);
    st2 = imfilter(st0,H,o.imfiltType);
    
    sshortind = st2<o.tth; % USE SMOOTHED TIMES FOR THRESHOLDING
    sh2(sshortind) = nan;
    st2(sshortind) = nan;
    
    shuffled_rs{reps} = sh2./st2;
end

%% DO THE AUTOCORRELATION (the time consuming part; 5 secs + per autocorrelation)

a = I2;
a = a*o.dsfs;
ac = corrcoef2(a,a,o.mincorrpixels);

shuffled_acs = cell(1,numreps);
for reps = 1:numreps
    sa = shuffled_rs{reps};
    sa = sa*o.dsfs;
    sac = corrcoef2(sa,sa,o.mincorrpixels);
    shuffled_acs{reps} = sac;
end

%%
opts = setdefaults([],'dozerpad',1,'nzer',2001,'detrend',0,'dotaper',1,'winType','hamming','taperStdDiv',8,'findPksInFFT',0);

x = R.xb2;
y = R.yb2;

dx = nanmean(diff(x));dy = nanmean(diff(y));
if dx~=dy
    error('use square pixels')
end

cx =ceil(size(ac,2)/2);cy =ceil(size(ac,1)/2);
[xx,yy] = meshgrid(-floor(size(ac,2)/2):floor(size(ac,2)/2),-floor(size(ac,2)/2):floor(size(ac,2)/2));
tapRadius = sqrt((xx.^2 + yy.^2));maxRad = max(tapRadius(:));
taper = coswineval(tapRadius,opts.winType,0,maxRad);

%%
b = ac;
b(isnan(b)) = 0;
b = b.*taper;

tmp = zeros(opts.nzer);
c =ceil(opts.nzer/2);
x = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dx;
y = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dy;
if isodd(size(b,1))
    tmp(c+[-floor(size(b,1)/2):floor(size(b,1)/2)],c+[-floor(size(b,2)/2):floor(size(b,2)/2)]) = b;
else
    tmp(c+[-floor(size(b,1)/2):floor(size(b,1)/2)-1],c+[-floor(size(b,2)/2):floor(size(b,2)/2)-1]) = b;
end
b = tmp;

spec = fftshift(fft2(b));
mag = abs(spec);
ps = 10* log10(abs((spec)).^2 );
%%

shuffled_mags = cell(1,numreps);
for reps = 1:numreps
    
    sb = shuffled_acs{reps} ;
    sb(isnan(sb)) = 0;
    sb = sb.*taper;
    
    stmp = zeros(opts.nzer);
    c =ceil(opts.nzer/2);
    x = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dx;
    y = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dy;
    if isodd(size(sb,1))
        stmp(c+[-floor(size(sb,1)/2):floor(size(sb,1)/2)],c+[-floor(size(sb,2)/2):floor(size(sb,2)/2)]) = sb;
    else
        stmp(c+[-floor(size(sb,1)/2):floor(size(sb,1)/2)-1],c+[-floor(size(sb,2)/2):floor(size(sb,2)/2)-1]) = sb;
    end
    sb = stmp;
    
    sspec = fftshift(fft2(sb));
    sspec = sspec(800:1200,800:1200);
    shuffled_mags{reps} = abs(sspec);
end

%%

avg_shuffled_mag = zeros(size(sspec,1),size(sspec,2),length(shuffled_mags));
for reps = 1:numreps;
    avg_shuffled_mag(:,:,reps) = shuffled_mags{reps}.^2;
end
avg_shuffled_mag = mean(avg_shuffled_mag,3)+std(avg_shuffled_mag,[],3);
%%
avg_shuffled_mag = zeros(size(sspec,1),size(sspec,2),length(shuffled_mags));
for reps = 1:numreps;
    avg_shuffled_mag(:,:,reps) = shuffled_mags{reps}.^2;
end
%%
noise_mag = zeros(size(sspec));
for rind = 1:size(sspec,1)
    for cind = 1:size(sspec,2)
        [~,~,CI] = ztest(avg_shuffled_mag(rind,cind,:),mean(avg_shuffled_mag(rind,cind,:)),...
            std(avg_shuffled_mag(rind,cind,:)),'alpha',0.001);
        noise_mag(rind,cind) = CI(2);
    end
end
avg_shuffled_mag2 = mean(avg_shuffled_mag,3);
%%
mag = mag(800:1200,800:1200);
cmag = mag.^2-avg_shuffled_mag2;
cmag(cmag < 0) = 0;
figure 
imagesc(cmag)
box off
title('mean noise removed')

cmag = mag.^2-noise_mag;
cmag(cmag < 0) = 0;
figure
imagesc(cmag)
box off
title('99.9% CI noise removed')


%%
fsx = 1/nanmean(diff(x));%samples per degree
fsy = 1/nanmean(diff(y));%samples per degree
max_freq_x = fsx/2;
max_freq_y = fsy/2;

xs = 0;ys = 0;if iseven(size(ps,2)), xs = -1;end;if iseven(size(ps,1)), ys = -1;end
x_freq_axis = [0:(floor(size(ps,2)/2))]/(floor(size(ps,2)/2))* max_freq_x; x_freq_axis = [-fliplr(x_freq_axis(2:end)) x_freq_axis(1:end+xs)];
y_freq_axis = [0:(floor(size(ps,1)/2))]/(floor(size(ps,1)/2))* max_freq_y; y_freq_axis = [-fliplr(y_freq_axis(2:end)) y_freq_axis(1:end+ys)];

x_periods = 1./x_freq_axis;
y_periods = 1./y_freq_axis;
xticksel = 1:10:length(x_freq_axis);
yticksel = 1:10:length(y_freq_axis);

dat = mag.^2;
%%
fig(101);
subplot 331;
imagesc(xb,yb,r);
xyl('X DVA','Y DVA',sprintf('Firing Rate Map'),14);
ax(1) = gca;

subplot 334;
imagescnan(x,y,b);axis xy;axis image;colormap fireprint;%freezeColors;colorbar
xyl('X DVA','Y DVA',sprintf('Tapered, Zero-Padded Autocorrelation'),14);
set(gca,'xlim',[-25 25],'ylim',[-25 25])
ax(4) = gca;

subplot 332;
fsel = nearest(x_freq_axis,-.3):nearest(x_freq_axis,0.3);
[xx,yy] = meshgrid(x_freq_axis,y_freq_axis);
imRads = sqrt(xx.^2+ yy.^2);
imagescnan(x_freq_axis(fsel),y_freq_axis(fsel),dat(fsel,fsel));
xyl('X Cycles/DVA','Y Cycles/DVA','2-D Power Spectrum',14)
axis square
ax(2) = gca;

subplot 333;
imagescnan(x_freq_axis(fsel),y_freq_axis(fsel),cmag(fsel,fsel));
xyl('X Cycles/DVA','Y Cycles/DVA','Cleaned Up 2-D Power Spectrum',14)
axis square
ax(3) = gca;

subplot 335
[fpcimg radii angles] = imgpolarcoord(dat(fsel,fsel),[],[0:0.5:359.5]/180*pi);radialSum = nansum(fpcimg,1);
xb3 = x_freq_axis(fsel);yb3 = fliplr(y_freq_axis(fsel));
dattmp = dat(fsel,fsel);
allradii = zeros(size(dattmp));allangle = allradii;
for yind = 1:size(dattmp,1)
    for xind = 1:size(dattmp,2)
        allradii(yind,xind) = norm([xb3(xind) yb3(yind)]);
        allangle(yind,xind) = 180/pi*atan2(yb3(yind),xb3(xind));
    end
end
radialSum = radialSum-min(radialSum);     radialSum = radialSum./max(radialSum);
ph = polar(angles,radialSum);set(ph,'color',BLUE)
xyl([],[],sprintf('2-D Power Spectrum\nin Polar Coordinates'),14)
ax(5) = gca;

subplot 338
plot(angles(nearest(angles,0):nearest(angles,pi))*180/pi,radialSum(nearest(angles,0):nearest(angles,pi)),'color',BLUE)
[pkInd pkValues ] = peakfinder([radialSum radialSum]);
peakAngles = angles(pkInd(pkInd<length(angles)))*180/pi ;
pkValues = pkValues(pkInd<length(angles));
set(gca,'xlim',[0 180],'xtick',rsig(peakAngles,0))
hold on;plot(peakAngles,pkValues,'x','color',RED,'markersize',10)
xyl('Grating Orientation (degrees)','Normalized Power','Periodic Band Power',14)
box off
ax(8) = gca;

subplot 336
[fpcimg radii angles] = imgpolarcoord(cmag(fsel,fsel),[],[0:0.5:359.5]/180*pi);radialSum = nansum(fpcimg,1);
xb3 = x_freq_axis(fsel);yb3 = fliplr(y_freq_axis(fsel));
dattmp = cmag(fsel,fsel);
allradii = zeros(size(dattmp));allangle = allradii;
for yind = 1:size(dattmp,1)
    for xind = 1:size(dattmp,2)
        allradii(yind,xind) = norm([xb3(xind) yb3(yind)]);
        allangle(yind,xind) = 180/pi*atan2(yb3(yind),xb3(xind));
    end
end
radialSum = radialSum-min(radialSum);     radialSum = radialSum./max(radialSum);
ph = polar(angles,radialSum);set(ph,'color',BLUE)
xyl([],[],sprintf('Cleaned Up 2-D Power Spectrum\nin Polar Coordinates'),14)
ax(6) = gca;

subplot 339
plot(angles(nearest(angles,0):nearest(angles,pi))*180/pi,radialSum(nearest(angles,0):nearest(angles,pi)),'color',BLUE)
[pkInd pkValues ] = peakfinder([radialSum radialSum]);
peakAngles = angles(pkInd(pkInd<length(angles)))*180/pi ;
pkValues = pkValues(pkInd<length(angles));
set(gca,'xlim',[0 180],'xtick',rsig(peakAngles,0))
hold on;plot(peakAngles,pkValues,'x','color',RED,'markersize',10)
xyl('Grating Orientation (degrees)','Normalized Power','Periodic Band Power',14)
box off
ax(9) = gca;


for atest = 1:length(peakAngles)
    [pkPower(atest) pkInd(atest)] = max(dattmp(allangle>(peakAngles(atest)-180-2) & allangle<(peakAngles(atest)-180+2)));
    pkFreqs(atest) = allradii(pkInd(atest));
end

pkFreqs = rsig(pkFreqs,3);
pkWavelength = 1./pkFreqs;
peakAngles = peakAngles;

relative_orientation = unique(diff(floor(peakAngles)));

subplot 337
bar(relative_orientation)
set(gca,'Xtick',[])
ylim([0 max(relative_orientation+10)])
box off
for pf = 1:length(pkFreqs)/2
    t= text(pf,relative_orientation(1)+10,[num2str(pkFreqs(pf)) 'cycles/dva']);
    set(t, 'HorizontalAlignment', 'center');
end
ylabel('Relative Orientation Angle (degrees)')

ax(7) = gca;

position_axes(gcf,ax,3,3,0,[2 2],10,[],[1 1])
