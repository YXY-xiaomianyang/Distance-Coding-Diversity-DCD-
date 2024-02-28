% Description:
%    This code is used to compute the distance codingdiversity (dcd).
% Parameter settings
s_win = 5;            % s_win is the window size
bin = 32;            % level is the number of gray levels

%%
% Output File Path
destination = 'E:\application_test';

%% image information
tic
[im,geo] = readgeoraster('E:\application_test\plant\test\test_plant.dat');
info = geotiffinfo('E:\application_test\plant\dcd_arrange\16_3_km_DCD_S.tif');
length_win = s_win*s_win;
a=floor(mean(im,3));
b=max(a,[],'all');
c=min(a,[],'all');
im = floor(a/((b-c)/bin)); 


%% DCD calculating

[n1,n2] = size(im);
DCD_img = zeros(n1,n2);

half = (s_win-1)/2;
jintival = half+1;
jendval = n2-half;

parfor i = half+1:n1-half
% tic
    for j = jintival:jendval
        S = DCD_S(im(i-half:i+half,j-half:j+half),length_win,bin);
        C = DCD_C(im(i-half:i+half,j-half:j+half),length_win,bin);
        SS = DCD_SS(im(i-half:i+half,j-half:j+half),length_win,bin);
        CS = DCD_CS(im(i-half:i+half,j-half:j+half),length_win,bin);
        CA = DCD_circle_5(im(i-half:i+half,j-half:j+half),length_win,bin);
        DCD_img(i,j) = mean([S C SS CS CA]);

    end

end

geotiffwrite([destination,'\',num2str(bin),'_',num2str(s_win),'_test_DCD.tif'],DCD_img,geo,'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

toc

%% DCD calculated by horizontal skip
function result = DCD_S(windw,length_win,bin)
if isempty(find(isnan(windw), 1))
    windw=windw.';
    Circle_win = [reshape(windw,1,[]),reshape(windw,1,[])];
    
    temp = [];

    for i = 1 : bin
        indx = find(Circle_win == i);
        if isempty(indx)
            continue
        end
        diff_indx = diff(indx);
        Disti = diff_indx(1:(length(diff_indx)+1)/2);
        tempi = count_hist(Disti);
        temp = [temp,tempi];
    end
    p = temp/length_win;
    result = sum(-p.*log(p));
else
    result = nan;
end
end

%% DCD calculated by vertical skip 
function result = DCD_C(windw,length_win,bin)
if isempty(find(isnan(windw), 1))
    
    Circle_win = [reshape(windw,1,[]),reshape(windw,1,[])];
    
    temp = [];

    for i = 1 : bin
        indx = find(Circle_win == i);
        if isempty(indx)
            continue
        end
        diff_indx = diff(indx);
        Disti = diff_indx(1:(length(diff_indx)+1)/2);
        tempi = count_hist(Disti);
        temp = [temp,tempi];
    end
    p = temp/length_win;
    result = sum(-p.*log(p));
else
    result = nan;
end
end

%% DCD calculated by horizontal serpentine
function result = DCD_SS(windw,length_win,bin)
[n,n]=size(windw);
if isempty(find(isnan(windw), 1))
    B=[];
    for i=1:n
        if mod(i,2)==1
            B=[B windw(i,:)];
        else
            B=[B fliplr(windw(i,:))];
        end
    end
    Circle_win = [B,B];
    
    temp = [];

    for i = 1 : bin
        indx = find(Circle_win == i);
        if isempty(indx)
            continue
        end
        diff_indx = diff(indx);
        Disti = diff_indx(1:(length(diff_indx)+1)/2);
        tempi = count_hist(Disti);
        temp = [temp,tempi];
    end
    p = temp/length_win;
    result = sum(-p.*log(p));
else
    result = nan;
end
end

%% DCD calculated by vertical serpentine
function result = DCD_CS(windw,length_win,bin)
windw=windw.';
[n,n]=size(windw);
if isempty(find(isnan(windw), 1))
    B=[];
    for i=1:n
        if mod(i,2)==1
            B=[B windw(i,:)];
        else
            B=[B fliplr(windw(i,:))];
        end
    end
    Circle_win = [B,B];
    
    temp = [];

    for i = 1 : bin
        indx = find(Circle_win == i);
        if isempty(indx)
            continue
        end
        diff_indx = diff(indx);
        Disti = diff_indx(1:(length(diff_indx)+1)/2);
        tempi = count_hist(Disti);
        temp = [temp,tempi];
    end
    p = temp/length_win;
    result = sum(-p.*log(p));
else
    result = nan;
end
end

%% DCD calculated by circle arrangement, take s_win as an example 
function result = DCD_circle_5(windw,length_win,bin)
if isempty(find(isnan(windw), 1))
    a=reshape(windw,1,[]);
    b=[a(13) a(18) a(23) a(24) a(19) a(25) a(20) a(14) a(15) a(10) a(9) a(5) a(4) a(8) a(3) a(2) a(7) a(1) a(6) a(12) a(11) a(16) a(17) a(21) a(22)];
    Circle_win = [b,b];
    
    temp = [];

    for i = 1 : bin
        indx = find(Circle_win == i);
        if isempty(indx)
            continue
        end
        diff_indx = diff(indx);
        Disti = diff_indx(1:(length(diff_indx)+1)/2);
        tempi = count_hist(Disti);
        temp = [temp,tempi];
    end
    p = temp/length_win;
    result = sum(-p.*log(p));
else
    result = nan;
end
end

%% Histogram Frequency Function
function c = count_hist(A)
sizeA = length(A);
temp=zeros(1,2500);
for pTa = 1:sizeA
    i = A(pTa)+1;
    temp(i) = temp(i)+1;
end
c = temp(temp>0);
end

