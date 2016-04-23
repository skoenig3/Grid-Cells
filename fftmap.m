function [fmap] = fftmap(R)
% Based on rmap: compute ratemap and autocorr input spike train, position values, and options by N Killian 110905

o = setdefaults(o,'getTAC',0,'useNaNfilt',0,'imfiltType',0,'useoldtimesforth',0,'mosersmth',0,'nbins',50,...
    'binningtype','specrangeExact','sqr',ceil(max(abs(Y))),'degperbin',0.5,'dogauss',1,'gsize',[13 13],'gsigma',3,...
    'spkth',0,'sth',1,'smtht',1,'timth',1,'tth',15,'dsfs',1e3,'mincorrpixels',15,'usexcorr2',0,...
    'keepSinglePixels',0,'Nrbin',0);

% 'symmetric', 'replicate','circular',0
rnanind = [];

if ~isfield(o,'sqr')&isfield(o,'squarerange'),o.sqr=o.squarerange;end
if o.mosersmth
    
    disp('using Moser method...')
    
    
    if isempty(strfind(o.smthmethod,'blackman'))
        Si = logical(S);  Sind = find(S);    Xi = X(Si);Yi = Y(Si);
    end
    
    switch o.smthmethod
        case 'blackman'
            %where S is a spike density/rate function for each position X,Y
            % this does spatial smoothing with a specified blackman window halfwidth
            % This is very similar to the NK Gaussian smoothing method using both time and spike smoothing
            halfwidth = o.blackhalf;%position values
            posi = 1;
            
            minval  = -o.sqr;maxval = o.sqr;
            nbins   = round([maxval-minval]/o.degperbin) + 1;
            dx = (maxval - minval) / (nbins-1);dy = dx;
            xr = minval:dx:maxval;   yr = minval:dy:maxval;
            
            lam = zeros(length(xr),length(yr));
            t2 = lam;
            
            wbar = waitbar(0,'calculating rate map...');  steps = length(xr)*length(yr);
            for xx = xr
                for yy = yr
                    %  weight = coswineval(dists,'blackman',0,halfwidth);
                    dists   = (((X-xx).^2 + (Y-yy).^2 ).^.5);
                    Stmp    = S(dists<halfwidth);
                    dists   = dists(dists<halfwidth);
                    dists   = dists/halfwidth/2 + .5;
                    weight  = (42 - 50*cos(2*pi*dists) + 8*cos(4*pi*dists))/100;%0 to 1
                    
                    lam(posi)  = sum(weight.*Stmp)/sum(weight);%weighted average spike count
                    t2(posi)   = sum(weight);
                    
                    waitbar(posi / steps)
                    posi = posi+1;
                end
            end
            close(wbar);
            
            xb = xr;yb = yr;
            I2 = lam;
            
            %             figure;imagesc(lam)
        case 'alldatabinned'
            %calculate the rate at each bin independently and if the time is too low then throw the bin out later, doesn't affect the other bins
            % this is the Moser method (Hafting 2005)
            disp('using alldatabinned method')
            posi = 1;
            minval  = -o.sqr;maxval = o.sqr;
            nbins   = round([maxval-minval]/o.degperbin) + 1;
            dx = (maxval - minval) / (nbins-1);dy = dx;
            xr = minval:dx:maxval;
            yr = minval:dy:maxval;
            
            lam = zeros(length(xr),length(yr));
            t2 = lam;
            
            wbar = waitbar(0,'calculating rate map...');  steps = length(xr)*length(yr);
            for xx = xr
                for yy = yr
                    smthspksum = sum(gausseval(( (Xi-xx).^2 + (Yi-yy).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));
                    smthtimsum = sum(gausseval(( (X-xx).^2 + (Y-yy).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));%leave as samples, convert to sec later
                    lam(posi)  = smthspksum/smthtimsum;
                    t2(posi) = smthtimsum;
                    waitbar(posi / steps)
                    posi = posi+1;
                end
            end
            close(wbar);
            
            xb = xr;yb = yr;
            I2 = lam;
            
        case 'alldatabinnednearest'
            closest = 5;
            posi    = 1;
            xr      = -o.sqr:o.degperbin:o.sqr;
            yr      = xr;
            lam     = zeros(length(xr),length(yr));
            t2 = lam;
            
            wbar = waitbar(0,'calculating rate map...');  steps = length(xr)*length(yr);
            for xx = xr
                for yy = yr
                    posdist = ((X-xx).^2 + (Y-yy).^2 ).^.5;spkdist = posdist(Si);
                    spkdist(posdist<=closest)=0;posdist = posdist(posdist<=closest);
                    if ~isempty(spkdist)&&~isempty(posdist)
                        smthspksum = nansum(gausseval(spkdist/o.smthfactor,0,o.mosgaussstd));
                        smthtimsum = nansum(gausseval(posdist/o.smthfactor,0,o.mosgaussstd));
                        lam(posi)  = smthspksum/smthtimsum;
                        t2(posi) = smthtimsum;
                    else
                        lam(posi) = 0;
                    end
                    waitbar(posi / steps)
                    posi = posi+1;
                end
            end
            close(wbar);
            
            xb = xr;yb = yr;
            I2 = lam;
            
        case 'justvisited' %estimate rate only at the points visited
            %has memory issues, and is not correct because rates are added at bins, not avg'd in hist3gc
            posi    = 1;
            xr      = sort(X); yr = sort(Y);
            lam     = zeros(length(xr),1);
            
            wbar    = waitbar(0,'calculating rate map...');steps = length(lam);
            for posi = 1:length(lam);
                smthspksum = nansum(gausseval(( (Xi-xr(posi)).^2 + (Yi-yr(posi)).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));
                smthtimsum = nansum(gausseval(( (X-xr(posi)).^2 + (Y-yr(posi)).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));
                lam(posi)  = smthspksum/smthtimsum;
                waitbar(posi / steps)
            end
            close(wbar)
            
            % now bin the values
            [h t r xb yb dx dy] = hist3gc(xr,yr,lam,o.nbins,o.nbins,o.binningtype,o.sqr,o.degperbin);
            I2 = r;
            t2 = t;
            
        case 'justvisitedloopbin' %estimate rate only at the points visited
            xr = X;yr = Y; lam = zeros(length(xr),1);
            
            minval  = -o.sqr;maxval = o.sqr;
            nbins   = round([maxval-minval]/o.degperbin) + 1;
            dx  = (maxval - minval) / (nbins-1);dy = dx;
            xrb = minval:dx:maxval;   yrb = minval:dy:maxval;
            
            lamb = cell(length(xrb),length(yrb));
            lam2 = zeros(length(xrb),length(yrb));
            time = zeros(length(xrb),length(yrb));
            
            wbar = waitbar(0,'calculating rate map...');steps = length(lam);
            for posi = 1:length(lam);
                xv = xr(posi);yv = yr(posi);
                positions(posi,:) = [xv yv];
                smthspksum = nansum(gausseval(( (Xi-xv).^2 + (Yi-yv).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));
                smthtimsum = nansum(gausseval(( (X-xv).^2 + (Y-yv).^2 ).^.5/o.smthfactor,0,o.mosgaussstd));
                lam(posi)  = smthspksum/smthtimsum;
                xbin = floor((xv-xrb(1)+dx/2)/dx)+1;
                ybin = floor((yv-yrb(1)+dy/2)/dy)+1;
                if xbin>0&ybin>0&xbin<nbins&ybin<nbins
                    lamb{xbin,ybin} = [lamb{xbin,ybin} lam(posi)];
                    time(xbin,ybin) = time(xbin,ybin)+1;
                end
                for xi = 1:numel(lam2)
                    lam2(xi) = nanmean(lamb{xi});
                end
                I2 = lam2;
                t2 = time;
                waitbar(posi / steps)
                
            end
            close(wbar);
            
    end
    
    if o.timth
        I2(t2<o.tth) = nan;
    end
    
    if o.thrate
        [pkrat pkratind] = max(I2(:));        rateth = o.thpercentage*pkrat;
        I2(I2<rateth) = 0;
        minarea = o.minarea/(o.degperbin^2);        % # bins^2
        stats = regionprops(I2>0,'Area','Centroid','PixelIdxList');
        for k = 1:length(stats)
            if stats(k).Area<=minarea
                I2(stats(k).PixelIdxList) = 0;
            end
        end
    end
    h0=[];r0 = [];t0=t2;x0=[];y0=[];z0=[];origind=[];h2 = [];
    
else
    
    %% do binning then smooth spikes and time, this is the Killian method
    % it's computationally much faster than the Moser-style methods,
    % but smooths over the kernel range after binning so is potentially less precise
    disp('using NK method...') 
    
    [h t r xb yb lb rb x0 y0 z0 origind] = hist3gc(X,Y,S,o.nbins,o.nbins,o.binningtype,o.sqr,o.degperbin);
   
    h0 = h;    t0 = t;%used for saving original data
    if o.spkth
        h(h<o.sth) = 0;
    end
    
    if ~o.smtht % don't smooth time
        t2=t+eps;
        r = h./t2;
        r(t==0)=nan; r0 = r;
        if o.timth
            r(t<o.tth) = 0;
            r0(t<o.tth) = nan;
        else
            rnanind = isnan(r);
        end
    else
        r0 = r;
    end
    
    r(rnanind) = 0;
    
    if o.dogauss
        H = fspecial('gaussian',o.gsize,o.gsigma);
        
        %  H = [0.0025 0.0125 0.0200 0.0125 0.0025;...
        %             0.0125 0.0625 0.1000 0.0625 0.0125;...
        %             0.0200 0.1000 0.1600 0.1000 0.0200;...
        %             0.0125 0.0625 0.1000 0.0625 0.0125;...
        %             0.0025 0.0125 0.0200 0.0125 0.0025;];%langston 2010 filter
       
        if ~o.smtht % Don't Smooth Time, Just Smooth The Binned Data
            if o.useNaNfilt
                disp('using nan filt');            I2 = nanfilt(r,H);
            else
                I2 = imfilter(r,H,o.imfiltType);
            end
            h2 = [];
            
            if o.timth
                shortind = t<o.tth;
                if o.keepSinglePixels
                    CC = bwconncomp(shortind,6);
                    stats = regionprops(CC,'Area','PixelIdxList');
                    for k = 1:length(stats)
                        if stats(k).Area<= 1 % TODO: MAKE GENERIC USING DEGPERBIN
                            shortind(stats(k).PixelIdxList) = 0;
                        end
                    end
                end
                %                 shortind is the same size as I2.
                I2(shortind) = nan;
            end
            
        else % Smooth Time and The Spikes Separately
            if o.useNaNfilt
                disp('using nan filt');   h2 = nanfilt(h,H);  t(t<eps) = nan;t2 = nanfilt(t,H);
            else
                h2 = imfilter(h,H,o.imfiltType); t2 = imfilter(t,H,o.imfiltType);
            end
            
            if o.timth
                if o.useoldtimesforth
                    shortind = t0<o.tth;
                else % USE SMOOTHED TIMES FOR THRESHOLDING
                    shortind = t2<o.tth;
                end
                
                % THIS WAS USED FOR THE PAPER BUT MAY NOT BE NECESSARY
                if o.keepSinglePixels
                    CC = bwconncomp(shortind,6);
                    stats = regionprops(CC,'Area','PixelIdxList');
                    for k = 1:length(stats)
                        if stats(k).Area<= 1 % TODO: MAKE GENERIC USING DEGPERBIN
                            shortind(stats(k).PixelIdxList) = 0;
                        end
                    end
                end
                h2(shortind) = nan;
                t2(shortind) = nan;
            end
            
            r = h2./t2;
            
            I2 = r;
        end
        
    else
        I2 = r;
    end
    
end
I2(rnanind) = nan;% restore the nans

% figure;imagescnan(I2);axis xy;axis image;pause
% THE RATEMAP IS NOW SAVED IN VARIABLE I2
% tic
if o.Nrbin>0.5
    disp('Digitizing the rate map')
    binLeft  = linspace(0,max(I2(:))+eps,o.Nrbin+1);
    binRight = binLeft(2:end);
    binLeft  = binLeft(1:end-1);
    binCenter = (binLeft+binRight)/2;
    for bini  = 1:numel(I2)
        if ~isnan(I2(bini))
%             I2(bini)
%             keyboard
            I2(bini) = binCenter(I2(bini)>=binLeft & I2(bini)<binRight);
        end
    end
end
% toc

%% DO THE AUTOCORRELATION

a = I2;
a = a*o.dsfs;

if ~isempty(o.bounds)
    a(:,xb<o.bounds(1)|xb>o.bounds(3)) = nan;
    a(yb<o.bounds(4)|yb>o.bounds(2),:) = nan;
end

% % % if isempty(strfind(o.binningtype,'normal'))
ac = corrcoef2(a,a,o.mincorrpixels);
if o.getTAC
    tac = corrcoef2(t2,t2,o.mincorrpixels);
else
    tac = [];
end
% % %     nanind = isnan(a);
ac2 = [];
% % %     a(nanind) = 0;% important,
% % %     ac2 = XCORR2x(a,a,'unbiased');
% % %     a(nanind) = nan;
% % % disp('sep. ac and ac2')
% % % else
% % %
% % %     ac2 = XCORR2x(a,a,'unbiased');
% % %
% % %     ac = ac2;
% % %     disp(' ac = ac2')
% % %
% % % end
% % % if o.usexcorr2
% % %     disp('using XCORR2x instead of corrcoef2()')
% % %     ac = ac2;
% % % end

xb2 = linspace(xb(1)*2,xb(end)*2,size(ac,2));
yb2 = linspace(yb(1)*2,yb(end)*2,size(ac,1));
nanind = isnan(ac);
ac(nanind) = 0;% important,

pks = [];acp= imregionalmax(ac);
[pks(:,2) pks(:,1)] = find(acp);
ac(nanind) = nan;

pkx = xb2(pks(:,1));
pky = yb2(pks(:,2));

R = var2field([],tac,nanind,h0,t0,r0,x0,y0,z0,origind,a,ac,ac2,pkx,pky,pks,xb,yb,xb2,yb2,t2,h2,0);

% % % figure;imagescnan(R.a);axis xy;axis image;pause


