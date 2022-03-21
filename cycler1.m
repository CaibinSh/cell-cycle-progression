function cell_cycle = cycler1(Hoechst,EdU,nuc_area,varargin)
%% CYCLER1(Hoechst,EdU,nuc_area,varargin) is designed to determine endpoint cell cycle phase
% INPUT:
%       Hoechst: a row vector representing integrated DNA content
%       measured by Hoechst staining
%       EdU: a row vector representing EdU signals
%       nuc_area: nuclear size (pixels)
%       varargin{1}(Optional): threshold for EdU signals, see EdgeDetecIF
%       varargin{2}(Optional): threshold for EdU signals, see EdgeDetecIF
% OUTPUT:
%       cell_cycle: 1--G1-phase; 2--S-phase; 3--G2-phase
%   
% reference: Gut, G., et al. Nat. Methods 12, 951?954 (2015).
% written by Caibin Sheng(shengcaibin@gmail.com), Loewer lab, TU Darmstadt

% check input
% set threshold for EdU signal determination, see EdgeDetecIF
if nargin == 3
    threshold = 0.003; 
    edge_range = 0.75;
elseif nargin == 4
    threshold = varargin{1};
    edge_range = 0.75;
elseif nargin == 5
    threshold = varargin{1};
    edge_range = varargin{2};
end
    
cycle = zeros(length(Hoechst),1);

% determine EdU positive cells (S-phase cells), see EdgeDetecIF
edge = EdgeDetecIF(EdU,threshold,edge_range);
S = find(EdU >= edge);G1_G2 = find(EdU < edge);
cycle(S)=2;

% determine G1- and G2-phase cells by unsupervised classification based on
% nuclear size and DNA content.
data_G1_G2 = [log10(Hoechst(G1_G2))',log10(nuc_area(G1_G2))'];
[idx,C] = kmeans(data_G1_G2,2,'Replicates',50);
G1 = idx==floor(sign(mean(data_G1_G2((idx==1),1))-mean(data_G1_G2((idx==2),1)))/2)+2;
    
cycle(G1_G2(G1))=1;
cycle(G1_G2(~G1))=3;

% data visualization -- scatter plot
Color = {[0 0.43 0.86],[0.57 0.29 0],[1 0.43 0.71]};
figure,
subplot(2,2,1)
for i=1:3
    scatter(log10(Hoechst(cycle==i)),log10(EdU(cycle==i)),8,...
            'filled','MarkerEdgeColor',Color{i},'MarkerFaceColor',Color{i},...
            'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.3)
    hold on
end

legend_and_partion = {['G1: ' num2str(round(sum(G1)/length(cycle)*100,3,'significant')),'%'],...
        ['S: ' num2str(round(length(S)/length(cycle)*100,3,'significant')),'%'],...
        ['G2: ' num2str(round(sum(~G1)/length(cycle)*100,3,'significant')),'%']}
xlim([5.4 6.4])
legend(gca,legend_and_partion,'Location','northwest')
legend('boxoff') 
ylabel('log10 (EdU signals(a.u.))')
xlabel('log10 (Hoechst intensity(a.u.))')
box on

% data visualization -- histogram
hold on
subplot(2,2,2),
for i =1:3
    h1 = histogram(Hoechst(cycle==i),min(Hoechst):(max(Hoechst)-min(Hoechst))/150:max(Hoechst));
    h1.FaceColor = Color{i};
    h1.EdgeColor = 'none';
    h1.FaceAlpha = 0.618
    xlim([min(Hoechst) max(Hoechst)])
    hold on
end

legend(gca,legend_and_partion)
ylabel('cell number')
xlabel('DNA content (Hoechst intensity (a.u.))')
cell_cycle = cycle;
end