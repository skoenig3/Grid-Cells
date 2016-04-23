% Compute Grid Cells
load('data.mat');

%%
spikes = find(data{2}.S == 1);
xpos = data{2}.X(spikes);
ypos = data{2}.Y(spikes);

plot(data{2}.X,data{2}.Y,'Color',[0.5 0.5 0.5])
axis tight
axis off
hold on

plot(xpos,ypos,'.r','markersize',12)
%%
xpos = round(xpos*24+400);
xpos(xpos > 800) = 800;
xpos(xpos < 1) = 1;
ypos = round(ypos*24+300);
ypos(ypos > 600) = 600;
ypos(ypos < 1) = 1;
ind = sub2ind([600 800],ypos,xpos);
grid = zeros(600,800);

for i = 1:length(ind);
    grid(ind(i)) = grid(ind(i))+1;
end
%%
filt = fspecial('gaussian',96,24);
gridfilt = imfilter(grid,filt);
%%
imagesc(log(abs(fftshift(fft2(gridfilt))).^2))

%%
xpos = round(data{2}.X*24+400);
xpos(xpos > 800) = 800;
xpos(xpos < 1) = 1;
ypos = round(data{2}.Y*24+300);
ypos(ypos > 600) = 600;
ypos(ypos < 1) = 1;
ind = sub2ind([600 800],ypos,xpos);
pos = zeros(600,800);
for i = 1:length(ind)
    pos(ind(i)) = pos(ind(i))+1;
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
%%
filt = fspecial('gaussian',8,3);
gridbinfilt = imfilter(gridbin,filt);
%%

grdbuffer = zeros(size(gridbinfilt));

template = [grdbuffer grdbuffer grdbuffer; 
            grdbuffer  gridbinfilt grdbuffer;
            grdbuffer grdbuffer grdbuffer];
        
shiftx = -size(gridbinfilt,2):size(gridbinfilt,2);
shifty = -size(gridbinfilt,1):size(gridbinfilt,1);
centerx = 99;
centery = 75;

indx = centerx-32:centerx+33;
indy = centery-24:centery+25;
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

