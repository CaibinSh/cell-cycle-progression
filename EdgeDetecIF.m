function edge = EdgeDetecIF(data,varargin)
%% EDGEDETECIF(data, varargin) is designed to determine positive signals for single cell immunofluorescence staining
%
% INPUT:
%       data: One-dimentional vector representing fluorescence siganls
%       varargin{1}(Optional): THRESHOLD, a threshold to determine edge
%       where signals increase at a high rate of spead
%       varargin{2}(Optional): EDGE_RANGE, the range where edge probably locates, e.g. 0.75 means
%edge on the right 3/4 of the data
% 
%OUTPUT:
%       edge: signal >= edge ---- positive; signal < edge ---- negative 
%   
% written by Caibin Sheng(shengcaibin@gmail.com), Loewer lab, MDC Berlin / TU Darmstadt

% check input
% default threshold and edge_range
    if nargin==1
        threshold = 0.003;
        edge_range = 0.75;
    elseif nargin==2
        threshold = varargin{1};
        edge_range = 0.75;
    elseif nargin==3
        threshold = varargin{1};
        edge_range = varargin{2};
    end
    
    % sort, smooth and normalize data
    data_sorted = sort(data);
    data_smooth = sgolayfilt(data_sorted,5,11);
    data_Znorm = zscore(data_smooth);
    
    % get increasing rate of signals
    FirDer = diff(data_Znorm);
    SecDer = diff(FirDer);
    
    % visualize increasing rate and threshold
    figure,
    subplot(2,1,1)
    plot(SecDer)
    a = find(SecDer>threshold);
    a(a<=(edge_range*length(data_sorted)))=[];
    title('increasing rate')
    ylim([-0.05 0.05])
    hold on, line([0.00 a(1)],[0.002 0.002])
    
    % visualize ranked cells and edge
    subplot(2,1,2)
    plot(data_sorted,'r.','MarkerSize',20)
    hold on,
    plot(data_smooth,'b');
    hold on,
    title('EdU signals')
    line([a(1) a(1)],[0 max(data_smooth)])
    xlabel('ranked single cells')
    ylabel('EdU (a.u.)')
    text(length(data)*0.618,max(data_smooth)*0.782,...
        {['EdU signals brighter than ',num2str(data_sorted(a(1)))];'considered as EdU positive cells'})
    edge = data_sorted(a(1));
end