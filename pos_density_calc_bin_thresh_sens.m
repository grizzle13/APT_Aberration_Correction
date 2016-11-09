function pos_density_calc_bin_thresh_sens

%% This function reads in position data from a standard IVAS .pos file and 
%  looks at density variation in precipitate volumes across the
%  reconstruction by sorting the ranged ions into user-defined bins,
%  discriminating by TOF, and fitting Gaussian curves to the resulting
%  distributions.
clear all
close all
clc

% Read in pos data
[posdata,pos_filename] = posread();

filename = pos_filename(1:end-4);

% Define run input parameters
[BinStartSize,BinEndSize,BinStepSize,zSkip,zTrunc,if_save] = binselect();

% Define peaks for TOF discrimination
% TOFselect utilizes a manual selection based on TOF spectrum
% staticTOF allows for specification of peak mass-to-charge limits
TOFposdata = TOFselect(posdata);
%TOFposdata = staticTOF(posdata,24.92,25.11,25.91,26.18,26.42,26.64);

csvoutput = zeros(length(BinStartSize:BinStepSize:BinEndSize)*3,13);

%Loop over defined bin sizes
n=1;
for BinSize = BinStartSize:BinStepSize:BinEndSize
    %Generate bins based on bin size
    [BinLimits,BinRanges,BinVectors,BinIndices] = makebins(posdata,...
        BinSize,zTrunc,zSkip);
            
    edges = cell(1,2);
    edges{1} = BinVectors{1};
    edges{2} = BinVectors{2};

    %Populate bins based on TOF of interest
    TOFcounts = posmap(TOFposdata,BinSize,BinLimits(5),BinIndices,edges);
    %Fit Gaussian and define minimum bin size to threshold
    nsigma = 2;
    for nsigma = 2:4
        figure(3*n-2)
        [cutoff,threshgof] = thresholdcounts(TOFcounts,nsigma);
        str = sprintf(['Thresholding for Bin Size = %1.2f, Threshold '...
            '= %u\\sigma'],BinSize,nsigma);
        title(str)
        if if_save == 'Y'
            ThresholdOut = sprintf('_Thresholding_bin%1.2f_%usigma',...
                BinSize,nsigma);
            ThresholdOut = strrep(ThresholdOut,'.','p');
            ThresholdOutName      = strcat(filename,ThresholdOut);
            print(ThresholdOutName,'-dpng');
        end
        %Repopulate bins based on threshold
        Threshcounts = posmap(TOFposdata,BinSize,BinLimits(5),...
            BinIndices,edges,cutoff);
        counts = posmap(posdata,BinSize,BinLimits(5),BinIndices,edges);
        figure(3*n-1)
        subplot(1,2,1)
        pcolor(edges{1},edges{2},counts(:,:,round(11/BinSize))'),...
            shading interp,axis square
        filtercounts = counts;
        filtercounts(Threshcounts==0)=0;
        subplot(1,2,2)
        pcolor(edges{1},edges{2},filtercounts(:,:,round(11/BinSize))'),...
            shading interp,axis square
        str = sprintf(['Density Hitmaps for Bin Size = %0.2f, Threshold'...
            '= %u\\sigma'],BinSize,nsigma);
        set(gcf,'NextPlot','add');
        axes;
        h = title(str);
        set(gca,'Visible','off');
        set(h,'Visible','on');
        if if_save == 'Y'
            HitmapOut = sprintf('_Hitmap_bin%1.2f_%usigma',BinSize,nsigma);
            HitmapOut = strrep(HitmapOut,'.','p');
            HitmapOutName      = strcat(filename,HitmapOut);
            print(HitmapOutName,'-dpng');
        end
        figure(3*n)
        str = sprintf(['Density Histograms for Bin Size = %0.2f, '...
            'Threshold = %u\\sigma'],BinSize,nsigma);
        [densification_factor,histogof] = histogen(counts,filtercounts);
        title(str)
        if if_save == 'Y'
            HistogramOut = sprintf('_Histogram_bin%1.2f_%usigma',...
                BinSize,nsigma);
            HistogramOut = strrep(HistogramOut,'.','p');
            HistogramOutName      = strcat(filename,HistogramOut);
            print(HistogramOutName,'-dpng');
        end
        csvoutput(n,1) = BinSize;
        csvoutput(n,2) = nsigma;
        csvoutput(n,3) = densification_factor;
        csvoutput(n,4) = threshgof.sse;
        csvoutput(n,5) = threshgof.rsquare;
        csvoutput(n,6) = threshgof.dfe;
        csvoutput(n,7) = threshgof.adjrsquare;
        csvoutput(n,8) = threshgof.rmse;
        csvoutput(n,9) = histogof.sse;
        csvoutput(n,10) = histogof.rsquare;
        csvoutput(n,11) = histogof.dfe;
        csvoutput(n,12) = histogof.adjrsquare;
        csvoutput(n,13) = histogof.rmse;
        n = n+1;
    end
    close all
end

if if_save == 'Y'
    csvwrite(strcat(filename,'_SensitivityAnalysis.csv'),csvoutput)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin custom functions

function [cutoff,gof] = thresholdcounts(TOFcounts,threshcoeff)
%%  Discriminates between matrix and precipitate volumes based on the bin
%   density distribution of solute ions based on the filtered position
%   data (TOFcounts) and a threshold density number of standard deviations
%   based on a number (threshcoeff) of standard deviations from the mean of
%   a Gaussian fit to the matrix bin density distribution.
    TOFpdfcount = TOFcounts(:);
    TOFpdfcount = TOFpdfcount(TOFpdfcount~=0);
    TOFpdfcount = sort(TOFpdfcount);
    TOFdist = histogram(TOFpdfcount,0:1:max(TOFpdfcount),...
        'normalization','pdf');
    TOFxdist = 0:1:length(TOFdist.Values)-1;
    TOFhisto = TOFdist.Values;

    % Currently uses Matlab curve-fitting add-on software package
    [xData, yData] = prepareCurveData( TOFxdist, TOFhisto );

    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    hold on
    gaussfit1 = fitresult.a1*exp(-((TOFxdist-fitresult.b1)...
        /fitresult.c1).^2);
    plot(TOFxdist,gaussfit1,'--r','LineWidth',1.5)
    
    cutoff=floor(fitresult.b1+threshcoeff*fitresult.c1);
    plot(cutoff,0,'og')
end

function [TOFcounts] = posmap(TOFposdata,BinSize,zmin,Nz,edges,varargin)
%%  Populate bins for with ions based input position data (TOFposdata) and
%   bin size parameters. varargin{1} assigns a density threshold, such that
%   bins with populations below this value are set to 0.
    for k = 1:Nz
        zslicedata = TOFposdata;
        zslicedata = zslicedata(:,zslicedata(3,:)>=(zmin+(k-1)*BinSize));
        zslicedata = zslicedata(:,zslicedata(3,:)<(zmin+k*BinSize));

        TOFcountmap = hist3([transpose(zslicedata(1,:)),...
            transpose(zslicedata(2,:))],'Edges',edges);
        if nargin > 5
            TOFcountmap(TOFcountmap<varargin{1})=0;
        end

        TOFcounts(:,:,k) = TOFcountmap;
    end
end

function [posdata,pos_filename] = posread()
%%  Reads in data from a standard .POS file generated by the IVAS software.
    disp('choose pos file');
    % prompt to read pos file
    [pos_filename,pos_path] = uigetfile({'*.pos'},'open');
    file_name=[pos_path pos_filename];

    pos_filename

    % reads position and mass-to-charge data from pos file
    fid = fopen(file_name, 'r');
    fseek(fid,0,-1);
    disp('loading position data');
    posdata=fread(fid, inf, 'float32', 0, 'b');
    posdata = reshape(posdata, [4 length(posdata)/4]); %Separates x, y, z, 
                                                 %m/n columns in .pos data
    disp('position data loaded');
end

%%%%%%%

function [BinStartSize,BinEndSize,BinStepSize,zSkip,zTrunc,if_save]...
    = binselect()
%%  Generates dialogue for inputting bin size parameters. Also allows for
%   truncation of the top and bottom parts of an APT data set to eliminate
%   some edge effects, and asks whether generated figures should be saved.
    disp('Input image resolution parameters (bin size)');
    prompt = {'start bin size (nm)','end bin size (nm)',...
        'bin step size (nm)','z depth to skip (nm)',...
        'z depth to truncate (nm)','save PNG + fig (Y/N)'};
    dlg_title = 'resolution parameters';
    num_lines = 1;
    def = {'.1','4','.05','8','16','N'};
    answer = (inputdlg(prompt,dlg_title,num_lines,def));

    BinStartSize = str2num(cell2mat(answer(1)));
    BinEndSize   = str2num(cell2mat(answer(2)));
    BinStepSize  = str2num(cell2mat(answer(3)));
    zSkip = str2num(cell2mat(answer(4)));
    zTrunc  = str2num(cell2mat(answer(5)));
    if_save = cell2mat(answer(6));

    if if_save(1)=='n'
        if_save='N'
    end
    if if_save(1)=='y'
        if_save='Y'
    end
end

%%%%%%%

function [TOFposdata] = TOFselect(posdata)
%%  Procedure for manually selecting solute ion peaks in the TOF spectrum
%   using a graphical user interface.
    figure
    histogram(posdata(4,:),0:0.01:ceil(max(posdata(4,:))))
    hold on
    xlabel('mass-to-charge (amu)')
    ylabel('counts')
    temp = input(['Adjust plot area to highlight desired ROI.'...
        ' Press enter when ready to continue']);
    disp('Click on the left boundary for the desired TOF range');
    title(['Click on the left boundary for the desired TOF range']);
    start_range=ginput(1)
    plot(start_range(1),0,'og')
    disp('Click on the right boundary for the desired TOF range');
    title(['Click on the right boundary for the desired TOF range']);
    end_range=ginput(1)
    plot(end_range(1),0,'og')
    if end_range<=start_range
        temp=start_range;
        start_range=end_range;
        end_range=temp;
        clear temp
    end
    TOFposdata = posdata;
    TOFposdata = TOFposdata(:,TOFposdata(4,:)>=start_range(1));
    TOFposdata = TOFposdata(:,TOFposdata(4,:)<=end_range(1));    

    choice = questdlg('Define additional ranges to include?', ...
        'Add ranges?','yes','no','no');

    switch choice
        case 'yes'
            flag_stop_ranges = 0;
        case 'no'
            flag_stop_ranges = 1;
    end

    while flag_stop_ranges<1
        figure(gcf)
        temp = input(['Adjust plot area to highlight desired ROI.'...
            ' Press enter when ready to continue']);

        disp('Click on the left boundary for the desired TOF range');
        title(['Click on the left boundary for the desired TOF range']);
        start_range=ginput(1)
        plot(start_range(1),0,'og')

        disp('Click on the right boundary for the desired TOF range');
        title(['Click on the right boundary for the desired TOF range']);
        end_range=ginput(1)
        plot(end_range(1),0,'og')

        if end_range<=start_range
            temp=start_range;
            start_range=end_range;
            end_range=temp;
            clear temp
        end

        TOFposdata1 = posdata;
        TOFposdata1 = TOFposdata1(:,TOFposdata1(4,:)>=start_range(1));
        TOFposdata1 = TOFposdata1(:,TOFposdata1(4,:)<=end_range(1));
        TOFposdata = [TOFposdata,TOFposdata1];
        clear TOFposdata1

        choice = questdlg('Define additional ranges to include?', ...
        'Add ranges?','yes','no','no');

        switch choice
            case 'yes'
                flag_stop_ranges = 0;
            case 'no'
                flag_stop_ranges = 1;
        end
    end    

    close(figure(gcf))

end

%%%%%%%

function [TOFposdata] = staticTOF(posdata,varargin)
%%  Allows for selection of solute ion peaks based on known peak positions
%   on the TOF spectrum (i.e., from a IVAS .rng file).
    k = (nargin-1)/2;
    for i = 1:k
        posrefdata = posdata;
        posrefdata = posrefdata(:,posrefdata(4,:)>=varargin{2*i-1});
        posrefdata = posrefdata(:,posrefdata(4,:)<=varargin{2*i});
        if i==1
           TOFposdata = posrefdata;
        else
           TOFposdata = horzcat(TOFposdata,posrefdata);
        end
    end
end

%%%%%%%

function [BinLimits,BinRanges,BinVectors,BinIndices] ...
    = makebins(posdata,BinSize,zTrunc,zSkip)
%%  Generates bin edges, limits and indices based on input parameters.
    xmin = min(posdata(1,:));
    xmax = max(posdata(1,:));
    ymin = min(posdata(2,:));
    ymax = max(posdata(2,:));
    %zmin = min(posdata(3,:));
    zmin = min(posdata(3,:))+zTrunc;
    %zmax = max(posdata(3,:));
    zmax = max(posdata(3,:))-zSkip;
    
    xrange = xmax - xmin;
    yrange = ymax - ymin;
    zrange = zmax - zmin;
    
    BinRanges = [xrange,yrange,zrange];
    
    xmin = xmin + rem(xrange,BinSize)/2;
    xmax = xmax - rem(xrange,BinSize)/2;
    ymin = ymin + rem(yrange,BinSize)/2;
    ymax = ymax - rem(yrange,BinSize)/2;
    zmin = zmin + rem(zrange,BinSize)/2;
    zmax = zmax - rem(zrange,BinSize)/2;
    
    BinLimits = [xmin,xmax,ymin,ymax,zmin,zmax];
    
    x = transpose(xmin:BinSize:xmax);
    y = transpose(ymin:BinSize:ymax);
    z = transpose(zmin:BinSize:zmax);
    
    BinVectors = cell(1,3);
    BinVectors{1} = x;
    BinVectors{2} = y;
    BinVectors{3} = z;
    
    Nx = length(x)-1;
    Ny = length(y)-1;
    Nz = length(z)-1;
    Ndat = length(posdata);
    
    BinIndices = [Nx,Ny,Nz,Ndat];
end

%%%%%%%

function [densification_factor,gof] = histogen(counts,filtercounts)
%%  Generates final histogram comparing bin densities in the matrix and in
%   the precipitate volumes and outputs the final densification factor.
    pdfcount = counts(:);
    pdfcount = pdfcount(pdfcount~=0);
    pdfcount = sort(pdfcount);

    filterpdfcount = filtercounts(:);
    filterpdfcount = filterpdfcount(filterpdfcount~=0);
    filterpdfcount = sort(filterpdfcount);

    hold on
    h_tot = histogram(pdfcount,0:1:max([filterpdfcount;pdfcount])+1);
    h_filt = histogram(filterpdfcount,...
        0:1:max([filterpdfcount;pdfcount])+1);
    h_net = h_tot.Values - h_filt.Values;
    xpdf = 0:1:length(h_net)-1;
    bar(xpdf+.5,h_net,1)
    h_filt = histogram(filterpdfcount,...
        0:1:max([filterpdfcount;pdfcount])+1);

    % Currently uses Matlab curve-fitting add-on software package
    [xData, yData] = prepareCurveData( xpdf, h_net );

    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    %opts.StartPoint = [9515 31 6.71367412863715];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    gaussfit2 = fitresult.a1*exp(-((xpdf-fitresult.b1)/fitresult.c1).^2);
    plot(xpdf,gaussfit2,'--g','LineWidth',1.5)

    bulk_dens = fitresult.b1;
    cluster_dens = mean(h_filt.Data);
    densification_factor = cluster_dens/bulk_dens
end