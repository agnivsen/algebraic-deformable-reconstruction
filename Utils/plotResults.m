function [rms] = plotResults(recons, gt)

    rms = RMSE3D(recons, gt);
    plot3(recons(:,1), recons(:,2), recons(:,3), 'o', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'b'); 
    hold on;
    plot3(gt(:,1), gt(:,2), gt(:,3), 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k'); 

    legend({'Reconstruction', 'G\textit{t}'}, 'Interpreter','latex', 'Location','northeast', 'FontSize', 8);
    title(['RMSE = ' num2str(rms)], 'Interpreter','latex');
    subtitle('Reconstructed vs. G\textit{t} pt.s', 'Interpreter','latex');
    axis on; grid on; axis equal;
    xlabel('X', 'Interpreter','latex');
    ylabel('Y', 'Interpreter','latex');
    zlabel('Z', 'Interpreter','latex');

end