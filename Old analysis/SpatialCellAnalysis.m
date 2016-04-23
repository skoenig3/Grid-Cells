% Compute Grid Cells
load('data.mat');
spikes = find(data{2}.S == 1);
xpos = data{2}.X(spikes);
ypos = data{2}.Y(spikes);

figure
plot(data{2}.X,data{2}.Y,'Color',[0.5 0.5 0.5])
hold on
plot(xpos,ypos,'.r','markersize',12)
hold off
axis tight
axis off
%%
xpos = round(xpos*24+400);
xpos(xpos > 800) = 800;
xpos(xpos < 1) = 1;
ypos = round(-ypos*24+300);
ypos(ypos > 600) = 600;
ypos(ypos < 1) = 1;
ind = sub2ind([600 800],ypos,xpos);
grid = zeros(600,800);

for i = 1:length(ind);
    grid(ind(i)) = grid(ind(i))+1;
end
%%
filt = fspecial('gaussian',128,24);
gridfilt = imfilter(grid,filt);
%%
time =zeros(600,800);
xspos = round(data{2}.X*24+400);
yspos = round(-data{2}.Y*24+300); 
xspos(xspos < 1) = 1; xspos(xspos > 800) = 800;
yspos(yspos < 1) = 1; yspos(yspos > 600) = 600; 
ind = sub2ind([600 800],yspos,xspos);
for i = 1:numel(time)
    time(i) = sum(ind == i); 
end
%%
binpos = bin2(pos,12,12,'lower','sum');
gridb = bin2(grid,12,12,'lower','sum');
gridbin = zeros(size(gridb));
for i = 1:size(binpos,1)
    for ii = 1:size(binpos,2);
        if binpos(i,ii) > 20
        gridbin(i,ii) = gridb(i,ii)/binpos(i,ii);
        end
    end
end
filt = fspecial('gaussian',8,3);
gridbinfilt = imfilter(gridbin,filt);
%%
[szr, szc] = size(gridbinfilt);
grdbuffer = zeros(szr, szc);

template = [grdbuffer grdbuffer grdbuffer; 
            grdbuffer  gridbinfilt grdbuffer;
            grdbuffer grdbuffer grdbuffer];
        
shiftx = -szc:szc;
shifty = -szr:szr;

if rem(size(template,2),2) == 0;
    centerx = size(template,2)/2;
else
    centerx = (size(template,2)+1)/2;
end
if rem(size(template,1),2) == 0;
    centery = size(template,1)/2;
else
    centery = (size(template,1)+1)/2;
end


indx = centerx-size(gridbinfilt,2)/2+1:centerx+size(gridbinfilt,2)/2;
indy = centery-size(gridbinfilt,1)/2+1:centery+size(gridbinfilt,1)/2;
xc = zeros(length(shifty),length(shiftx));
for i = 1:length(shiftx);
    for ii = 1:length(shifty);
        c = corrcoef(gridbinfilt,template(indy+shifty(ii),indx+shiftx(i)));
        xc(ii,i) = c(2);
    end
end
%%
binpos = bin2(pos,12,12,'lower','sum');
shuffledpower = zeros(size(gridbinfilt));
filt = fspecial('gaussian',16,3);
for i = 1:100
    sgrid = grid(randperm(numel(grid)));
    sgrid = reshape(sgrid,size(grid));
    binned = bin2(sgrid,12,12,'lower','sum');
    binshuff = zeros(size(gridb));
    for i = 1:size(binpos,1)
        for ii = 1:size(binpos,2);
            if binpos(i,ii) > 20
                binshuff(i,ii) = gridb(i,ii)/binpos(i,ii);
            end
        end
    end

    binshuff = imfilter(binshuff,filt);
    power = fftshift(fft(binshuff));
    shuffledpower = shuffledpower + power;
end
%%
gridfftpower = abs(fftshift(fft2(gridbinfilt)));
imagesc(gridfftpower-shuffledpower)

