function Heatmap(u, U)
    figure
    imagesc([-5, 5], [-5, 5], u)
    colorbar
    axis square
    colormap(othercolor('OrRd9'))
%     colormap('winter')
end