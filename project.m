%% Author: Nathaniel Mathews
%% Initialize

clear all; close all;
% Put this file in a folder named 'stats project', and put
% the directory to that folder here, eg
DIRECTORY = 'D:/Documents/Matlab/stats project';
addpath('stats project');

b0 = @(x) 7*sin(0.0172*x-2.68);

FIGURE = 1;

%% Set years (test bed)
years = {'2007','2008','2009'};

%% Set years (validation)
years = {'1999','2000'};


%% Read in flare data
flare_dates = [0];
n = 0;
for y = 1:size(years,2)
    f = dir(sprintf('%s/%s/flare-events_*.txt',DIRECTORY,years{y}));
    for j=1:size(f)
        % flare-event files formatted as '5s%6s%6s%5s%5s%14s%10s%s11%13s%5s%5s%4s%6s'
        % goes-xrs and flare-events both have the date 6:11, patrol as 7:12
        fid = fopen(sprintf('%s/%s/%s',DIRECTORY,years{y},f(j).name));
        tline = fgetl(fid);
        tline = fgetl(fid); % Use only for h-alpha flare tables
        while size(tline,2) >= 2
            n = n + 1;
            flare_dates(n) = datenum(tline(6:11),'yymmdd');
            tline = fgetl(fid);
            tline = fgetl(fid); % Use only for h-alpha flare tables
        end
        fclose(fid);
    end
end

clear fid n tline j f;

%% Read  in sunspot data
sunspots = [0,0,0];
n = 0;
for y = 1:size(years,2)
    f = dir(sprintf('%s/%s/usaf_solar-region-reports_*.txt',DIRECTORY,years{y}));
    for j=1:size(f)
        % details available at http://www.ngdc.noaa.gov/stp/space-weather/
        % solar-data/solar-features/sunspot-regions/usaf_mwl/documentation/
        % usaf_sunspot-region-reports_format.txt
        
        fid = fopen(sprintf('%s/%s/%s',DIRECTORY,years{y},f(j).name));
        tline = fgetl(fid);
        while size(tline,2) >= 2
            if ~any(sunspots(:,1) == datenum(tline(3:8),'yymmdd')) && ...
                    sum(isstrprop(tline(49:52),'digit')) && ...
                    sum(isstrprop(tline(44:45),'digit'))
                n = n + 1;
                sunspots(n,1) = datenum(tline(3:8),'yymmdd');
                sunspots(n,2) = str2double(tline(44:45));
                sunspots(n,3) = str2double(tline(49:52))*0.000001;
            end
            tline = fgetl(fid);
        end
        fclose(fid);
    end
end

clear f fid j n tline;

%% Read in Mount Wilson Magnetograms

resolution = 200;
N = resolution / 10;
coronal_dates = [0];
variances = [];
means     = [];
modes     = [];
overall_index = 0;
for y = 1:size(years,2)
    f = dir(sprintf('%s/%s/mwil_vmgcv_fp_*',DIRECTORY,years{y}));
    w = waitbar(0, sprintf('reading in %s data',years{y}));
    for coronograph = 1:size(f)
        overall_index = overall_index+1;
        waitbar(coronograph/size(f,1))
        if size(f(coronograph).name,2) == 24
            coronal_dates(overall_index) = datenum(f(coronograph).name(15:20),'yymmdd');
        elseif size(f(coronograph).name,2) == 26
            coronal_dates(overall_index) = datenum(f(coronograph).name(17:22),'yymmdd');
        else
            disp(sprintf('Unexpected file name length, name = "%s"', ...
                f(coronograph).name));
        end
        cdata = imread(sprintf('%s/%s/%s',DIRECTORY,years{y},f(coronograph).name));
        cdata = squeeze(cdata(:,:,1)); % necessary for a few wonky jpgs
        resize = resolution/min(size(cdata));
        cdata = imresize(255/max(max(cdata))*cdata,resize);

        snips = zeros(floor(size(cdata,1)/N)*floor(size(cdata,2)/N),N,N);
        j = 1; n = 1; while j < size(cdata,1)-N
            k = 1;
            while k < size(cdata,2)-N
                snips(n,:,:) = cdata(j:j+N-1,k:k+N-1);
                k = k+N;
                n = n+1;
            end
            j = j+N;
        end
        j = 1; while j <= size(snips,1)
            v = reshape(squeeze(snips(j,:,:)),1,[]);
            if (sum(v >= max(max(cdata))-2) >= size(v,2)/4) || ...
                (sum(v <= min(min(cdata))+2) >= size(v,2)/4)
                snips(j,:,:) = [];
            else
    %             if coronograph == 10
    %                 figure(FIGURE); FIGURE = FIGURE + 1;
    %                 imagesc(squeeze(snips(j,:,:)))
    %                 colorbar;
    %             end
                
                means(overall_index, j) = mean(double(v) ...
                    *15/256-7.5+b0(coronal_dates(overall_index)-datenum('01/01/99')));
                modes(overall_index, j) = mode(v ...
                    *15/256-7.5+b0(coronal_dates(overall_index)-datenum('01/01/99')));
                variances(overall_index, j) = var(v ...
                    *15/256-7.5+b0(coronal_dates(overall_index)-datenum('01/01/99')));
                j = j+1;

            end
        end
    end
    close(w)
    
end
sortedVars = sort(variances,2,'descend');
sortedMeans = sort(means,2,'descend');
sortedModes = sort(modes,2,'descend');
colMax = 41;

clear cdata coronograph f j k n overall_index resize snips v w y;
clear sortedVars sortedMeans sortedModes;

%% Construct X Matrix

flare_counts = [0];
r_sunspots  = [0,0];
r_vars      = zeros(1,colMax);
r_means     = zeros(1,colMax);
r_years     = [0];
n = 0;
% for y = max(min(coronal_dates),min(sunspots(:,1))): ...
%         min(max(coronal_dates),max(sunspots(:,1)))
%     if any(coronal_dates == y) && any(sunspots(:,1) == y)
%         n = n+1;
%         flare_counts(n) = sum(flare_dates == y);
%         r_sunspots(n,1) = sunspots(sunspots(:,1) == y,2);
%         r_sunspots(n,2) = sunspots(sunspots(:,1) == y,3);
%         r_vars(n,:) = variances(coronal_dates == y,1:colMax);
%         r_means(n,:) = abs(means(coronal_dates == y,1:colMax));
%         r_years(n) = y;
%     end
% end
n = 0;
for y = coronal_dates
    n = n + 1;
    flare_counts(n) = sum(flare_dates == y);
    if any(sunspots(:,1) == y)
        r_sunspots(n,1) = sunspots(sunspots(:,1) == y,2);
        r_sunspots(n,2) = sunspots(sunspots(:,1) == y,3);
    else
        r_sunspots(n,1) = 0;
        r_sunspots(n,2) = 0;
    end
    r_vars(n,:) = variances(coronal_dates == y,1:colMax);
    r_means(n,:) = abs(means(coronal_dates == y,1:colMax));
    r_years(n) = y;
end
diff_vars = [zeros(1,colMax);r_vars(2:end,:)-r_vars(1:end-1,:)];
[~,sorting] = sort(r_means,2,'descend');
sort_vars = zeros(size(r_vars));
sort_means = zeros(size(r_means));
sd_vars = zeros(size(r_vars));
for j = 1:size(r_means,1)
    sort_vars(j,:)  = r_vars(   j,sorting(j,:));
    sort_means(j,:) = r_means(  j,sorting(j,:));
    sd_vars(j,:)    = diff_vars(j,sorting(j,:));
end

% disp('Basic_glm')
% X = [r_sunspots,mean(r_vars,2),mean(r_means,2)];
% 
% disp('Unsorted_glm')
% X = [r_sunspots,r_vars,r_means];
% 
disp('Sorted_glm')
X = [r_sunspots,sort_vars,sort_means];
% 
% disp('Memory_glm')
% X = [r_sunspots,sort_vars,sort_means,sd_vars];

%% Fit glm
% 
disp('Apply Principal Components')
[xcoeff,xscore,xlatent,~,xexplained] = pca( X );
X = X*xcoeff;
X = X(:,1:85);
% 
% figure(FIGURE); FIGURE = FIGURE+1;
% semilogy(xlatent)
% xlabel('Principal component')
% ylabel('Singular Value')

glm = fitglm(X,flare_counts, 'Distribution','poisson');

% figure(FIGURE); FIGURE = FIGURE+1;
% plotResiduals(glm,'fitted')
% 
% figure(FIGURE); FIGURE = FIGURE+1;
% hold on;
% plot(flare_counts, flare_counts);
% plot(flare_counts, glm.Fitted.Response, 'x');
% legend('True values','Fitted values')
% xlabel('Flare counts')
% ylabel('Flare counts')
% hold off;
% 
% figure(FIGURE); FIGURE = FIGURE+1;
% plotDiagnostics(glm,'contour')
% 
% figure(FIGURE); FIGURE = FIGURE+1;
% plotDiagnostics(glm,'cookd')
% 
% figure(FIGURE); FIGURE = FIGURE+1;
% plotResiduals(glm,'lagged')


%% Validation

% X = X*xcoeff; X = X(:,1:85);

[predictions, ci] = predict(glm3,X);

p = mean(flare_counts' < ci(:,2) & flare_counts' > ci(:,1))

pci  = [predictions - ci(:,1), ci(:,2) - predictions];

figure(FIGURE); FIGURE = FIGURE+1;
hold on;
title('1999')
plot(r_years(1:22),flare_counts(1:22))
errorbar(r_years(1:22),predictions(1:22),pci(1:22,1),pci(1:22,2),'x')
legend('True flare counts','predicted flare counts')
xlabel('Date')
ylabel('Flare counts')
axis([-inf inf 0 50])
hold off;

figure(FIGURE); FIGURE = FIGURE+1;
hold on;
title('2000')
plot(r_years(23:end),flare_counts(23:end))
errorbar(r_years(23:end),predictions(23:end),pci(23:end,1),pci(23:end,2),'x')
legend('True flare counts','predicted flare counts')
xlabel('Date')
ylabel('Flare counts')
axis([-inf inf 0 25])
hold off;

