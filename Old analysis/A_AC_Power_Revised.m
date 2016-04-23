clar
% load respcell;
% % % load validunits;respcell = strcat(validunits(:,1),{' '},validunits(:,2));
% load manualsel;
dir('C:\Users\seth.koenig\Documents\MATLAB\Grid Cells');genfigs = 1;

% gauss, hamming, hann, blackman, flattopwin
% clear w;figure;
% w(:,1) = coswineval(1:51,'hann',26,25);
% w(:,2) = coswineval(1:51,'hamming',26,25);
% w(:,3) = coswineval(1:51,'blackman',26,25);
% w(:,4) = coswineval(1:51,'flattopwin',26,25);
% plot(w)
% legend('hann','hamming','blackman','flattopwin')
% set(gca,'xlim',[1 51],'ylim',[-.1 1]);box off;

opts = setdefaults([],'dozerpad',1,'nzer',2001,'detrend',0,'dotaper',1,'winType','hamming','taperStdDiv',8,'findPksInFFT',0);

% for k = 1:size(respcell,1)
for k = 5
    %     load(['R:\Buffalo Lab\Nathan\backup\RateMaps\' respcell{k,1} 'ratemapFinal2.mat'])
         load('R:\Buffalo Lab\Nathan\backup\RateMaps\MP111213 sig007cratemapFinal2.mat')
      %         %load('R:\Buffalo Lab\Nathan\backup\RateMaps\MP111213 sig004aratemapFinal2.mat')
     %load('R:\Buffalo Lab\Nathan\backup\RateMaps\TT091218 sig003aratemapFinal2.mat')
    
    % PS has same dimensions as the original image, but just take half.
    R.ac(isnan(R.ac)) = 0;b = R.ac; b0 = b;x = R.xb2; y = R.yb2;
    % R.a(isnan(R.a)) = 0;b = R.a; b0 = b;x = R.xb2; y = R.yb2;
    
    %     wDeg = 24;  %size of image (in degrees)
    %     nPix = 95;  %resolution of image (pixels);
    %     [x,y] = meshgrid(linspace(-wDeg/2,wDeg/2,nPix+1));
    %     x = x(1:end-1,1:end-1);y = y(1:end-1,1:end-1);
    %     %deg (counter-clockwise from vertical)
    %     sf = 0.15; %spatial frequency (cycles/deg)
    %     orientation = -15;orientation2 = 50;orientation3 = 80;
    % %     orientation = -30;orientation2 = 30;orientation3 = 90;
    %     ramp = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;
    %     ramp2 = cos(orientation2*pi/180)*x - sin(orientation2*pi/180)*y;
    %     ramp3 = cos(orientation3*pi/180)*x - sin(orientation3*pi/180)*y;
    %     % grating = sin(2*pi*sf*ramp).*sin(2*pi*sf*ramp2).*sin(2*pi*sf*ramp3);
    %     grating = sin(2*pi*sf*ramp) + sin(2*pi*sf*ramp2) + sin(2*pi*sf*ramp3);
    %     x = linspace(-wDeg/2,wDeg/2,nPix);y = x;
    %     b = grating;
    %
    if size(b,2)~=size(b,1),error('not everything is compatible with non-square matrices yet!');end
    
    if opts.dotaper
        
        if strcmp(opts.winType,'gauss')
            taper = fspecial('gauss',size(b),size(b,2)/opts.taperStdDiv); taper = taper/max(taper(:));
        else
            cx =ceil(size(b,2)/2);cy =ceil(size(b,1)/2);[xx yy] = meshgrid(-floor(size(b,2)/2):floor(size(b,2)/2),-floor(size(b,2)/2):floor(size(b,2)/2));
            tapRadius = sqrt((xx.^2 + yy.^2));maxRad = max(tapRadius(:));
            taper = coswineval(tapRadius,opts.winType,0,maxRad);
        end
        %         figure;imagesc(taper);        pause
        
        b = b.*taper;
        
    end
    if opts.detrend
        b = b-nanmean(b(:));
    end
    
    dx = nanmean(diff(x));dy = nanmean(diff(y));
    if dx~=dy, error('use square pixels');end
    
    if opts.dozerpad
        tmp = zeros(opts.nzer);c =ceil(opts.nzer/2);x = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dx;y = [-floor(opts.nzer/2):floor(opts.nzer/2)] * dy;
        if isodd(size(b,1))
            tmp(c+[-floor(size(b,1)/2):floor(size(b,1)/2)],c+[-floor(size(b,2)/2):floor(size(b,2)/2)]) = b;
        else
            tmp(c+[-floor(size(b,1)/2):floor(size(b,1)/2)-1],c+[-floor(size(b,2)/2):floor(size(b,2)/2)-1]) = b;
        end
        b = tmp;
    end
    
    
    spec = fftshift(fft2(b));
    pha = angle(spec)/pi*180;
    mag = abs(spec);
    ps = 10* log10(abs((spec)).^2 );
    
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
    
    
    %%
    ax = [];
    fig(101);
    %         title(sprintf('%g: %s',k,respcell{k}),'fontsize',16)
    
    subplot 221;
    imagescnan(x,y,b);axis xy;axis image;colormap fireprint;%freezeColors;colorbar
    xyl('X DVA','Y DVA',sprintf('Tapered, Zero-Padded Autocorrelation'),14);
    set(gca,'xlim',[-25 25],'ylim',[-25 25])
%     set(gca,'xlim',[-33 33],'ylim',[-33 33])
%         set(gca,'xlim',[-15 15],'ylim',[-15 15])

    dat = mag.^2;
    %         dat = mag;
    ax(1) = gca;
    
    subplot 222;
    % % %     ax(2) = gca;cla(ax(2));
    
    %     fsel = 1:size(mag,2);
    %     fsel = nearest(x_freq_axis,-.5):nearest(x_freq_axis,0.5);
    %     fsel = nearest(x_freq_axis,-.4):nearest(x_freq_axis,0.4);
    fsel = nearest(x_freq_axis,-.3):nearest(x_freq_axis,0.3);
    %     fsel = nearest(x_freq_axis,-.2):nearest(x_freq_axis,0.2);
    
    [xx yy] = meshgrid(x_freq_axis,y_freq_axis);
    imRads = sqrt(xx.^2+ yy.^2);
    
    %     tooWideSel = nearest(x_freq_axis,-.03):nearest(x_freq_axis,0.03);inhib = 1./fspecial('gauss',[length(tooWideSel) length(tooWideSel)],length(tooWideSel)/.1);inhib = inhib-min(inhib(:));inhib = inhib/max(inhib(:));%         dat(tooWideSel,tooWideSel) = dat(tooWideSel,tooWideSel).*inhib;
    
    %     tooWideSel = imRads<0.03;dat(tooWideSel) = 0;
    
    
    imagescnan(x_freq_axis(fsel),y_freq_axis(fsel),dat(fsel,fsel));
    %         imagescnan(x_freq_axis(fsel),y_freq_axis(fsel),pha(fsel,fsel));
    
%     colormap fireprint;    colorbar; hold on;    axis image;
    xyl('X Cycles/DVA','Y Cycles/DVA','2-D Power Spectrum',14)
%     set(gca,'xtick',x_freq_axis(xticksel),'ytick',y_freq_axis(yticksel),'xticklabel',rsig(x_periods(xticksel),2),'yticklabel',rsig(y_periods(yticksel),2))
    
    
    if opts.findPksInFFT
        % pks = FastPeakFind(dat(fsel,fsel),20,fspecial('gaussian', 5 , 1), 2, 0);%7c
        pks = FastPeakFind(dat(fsel,fsel),10,[], 2, 0);%peaks in the 2D fft is a little futile, only works with perfect data
        for pkind = 1:size(pks,1)
            plot(x_freq_axis(fsel(pks(pkind,2))),y_freq_axis(fsel(pks(pkind,1))),'bx','markersize',3)
            plot(x_freq_axis(fsel(pks(pkind,2))),y_freq_axis(fsel(pks(pkind,1))),'ro','markersize',5)
        end
    end
    ax(2) = gca;
    
    
    subplot 223
    %     [fpcimg radii angles] = imgpolarcoord(dat(fsel,fsel),[],[0:0.5:359.5]/180*pi);
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
    % % %     [Y ,angles, radialSum] = binanalslide1D(allangle(:)',[-180:1:180],1,1,dattmp(:)'); angles = angles/180*pi;
    
    radialSum = radialSum-min(radialSum);     radialSum = radialSum./max(radialSum);
    
    ph = polar(angles,radialSum);set(ph,'color',BLUE)
    xyl([],[],sprintf('2-D Power Spectrum\nin Polar Coordinates'),14)
    
    ax(3) = gca;
    
    
    
    
    subplot 224
    plot(angles(nearest(angles,0):nearest(angles,pi))*180/pi,radialSum(nearest(angles,0):nearest(angles,pi)),'color',BLUE)
    [pkInd pkValues ] = peakfinder([radialSum radialSum]);
    peakAngles = angles(pkInd(pkInd<length(angles)))*180/pi ;
    pkValues = pkValues(pkInd<length(angles));
    set(gca,'xlim',[0 180],'xtick',rsig(peakAngles,0))
    hold on;plot(peakAngles,pkValues,'x','color',RED,'markersize',10)
    xyl('Grating Orientation (degrees)','Normalized Power','Periodic Band Power',14)
    
    for atest = 1:length(peakAngles)
        [pkPower(atest) pkInd] = max(dattmp(allangle>(peakAngles(atest)-180-2) & allangle<(peakAngles(atest)-180+2) & allradii>0.03 & allradii<0.3));
        [pkFreqs(atest,2) pkFreqs(atest,1)] = ind2sub(size(allangle),pkInd);
        pkFreqs(atest,2) =  yb3(pkFreqs(atest,2));
        pkFreqs(atest,1) =  xb3(pkFreqs(atest,1));
        
    end
    %     pkFreqs =   %FIXME
    pkWavelength = 1./pkFreqs;
    peakAngles = peakAngles;
    
    wrapped = wrapToPi(peakAngles*pi/180);wrapped(wrapped<0) = wrapped(wrapped<0)+pi;
    wrapped = unique(rsig(wrapped,3));
    
    ax(4) = gca;
    
    
    %         pause
    if genfigs
        %         dispfig;
        position_axes(gcf,ax,2,2,0,[2 2],12,[],[1.5 1])
%         export_fig([figdir sprintf('%s_%s',leadz(k,2),respcell{k}) '.png'])
%         export_fig([figdir sprintf('%s_%s',leadz(k,2),respcell{k}) '.pdf'])
    end
    
    results(k) = var2field([],peakAngles,pkValues,pkWavelength,pkFreqs,wrapped,angles,radialSum,0);
    
end

% save2(['D:\Dropbox\fft2dat_' opts.winType '.mat'],results,respcell,opts)



