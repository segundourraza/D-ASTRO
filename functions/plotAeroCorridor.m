%% FLIGHT CORRIDOR PLOT

%Parameters for making plots
Markers = {'.','+','x','*','o','v','d','^','s','>','<'};
Line_format={'--','-',':','-.'};
Color = {'b','r'};
Counter=0;
Markermax=length(Markers);
Line_formatmax=length(Line_format);
Colormax=length(Color);

density_label = {'normal', 'low', 'high'};
BC_range = simInputs.BC_range;
nBC = length(BC_range);

if length(BC_range)>1
    hold on
    if simInputs.completeCorridor
        lb=Gammas(1:nBC, [1,3,5])'*180/pi; % lower bound flight path angles stored here
        ub=Gammas(1:nBC, [2,4,6])'*180/pi; % upper bound flight path angles stored here
        bisection = (lb(3,:)+ub(2,:))/2;
        patch([BC_range, fliplr(BC_range)], [lb(1,:), fliplr(ub(1,:))], 'y','FaceAlpha',.7)
        index = find(lb(3,:) < ub(2,:));
        patch([BC_range(index) fliplr(BC_range(index))], [ub(2,index) fliplr(lb(3,index))], 'g','FaceAlpha',.7)

        for i = 1 : 3
            plot(BC_range,lb(i,:),strcat(Line_format{mod(i,Line_formatmax)+1},Color{mod(1,Colormax)+1}),'LineWidth',2)
            plot(BC_range,ub(i,:),strcat(Line_format{mod(i,Line_formatmax)+1},Color{mod(2,Colormax)+1}),'LineWidth',2)
        end
        legLabel = {'','','normal: minimum','normal: maximum', 'low: minimum','low: maximum','high: minimum','high: maximum'};
    else
        lb=Gammas(1:nBC, 1)'*180/pi; % lower bound flight path angles stored here
        ub=Gammas(1:nBC, 2)'*180/pi; % upper bound flight path angles stored here    
        bisection = (lb+ub)/2;
        patch([BC_range fliplr(BC_range)], [ub fliplr(lb)], 'g','FaceAlpha',.7, 'DisplayName', 'Robust Corridor')
        legLabel = {};
    end
    maxGamma = max(ub,[],"all")*0.95;
    minGamma = min(lb, [], 'all')*1.05;

    plot(BC_range, bisection, '--k', 'LineWidth',1, DisplayName='Corrdior Bisector')
    str = {'Surface','Collision'};
    text(BC_range(1)*1.2,minGamma*0.95,str,'FontSize',18)
    str2 = {'Unbounded','Orbit'};
    text(BC_range(end)*0.5,maxGamma*1.1,str2,'FontSize',18)
    hold off
    
    grid minor
    ylabel("\gamma (\circ)")
    xlabel('\beta')
    title('Aerocapture Flight Corridor')

    xlim([BC_range(1), BC_range(end)])
    ylim([minGamma, maxGamma])
    set(gca,'GridAlpha',0.4,'MinorGridAlpha',0.7);
    set(gca, 'FontSize', 18)
    leg1 = legend(legLabel,'Location','eastoutside', 'fontsize', 12);
end