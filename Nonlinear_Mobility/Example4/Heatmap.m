function Heatmap(u)
    figure
    imagesc(u, [-1, 1])
    colorbar
    axis square
    colormap(othercolor('GnBu9'))
%     colormap('winter')
end