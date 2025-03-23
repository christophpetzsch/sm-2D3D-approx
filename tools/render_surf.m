function render_surf(X)
% This code is adopted from Marvin Eisenberg.

V = X.VERT;
F = X.TRIV;


%% Plot mesh
t = tsurf(F,V, 'FaceColor', '#DCDDD8', 'FaceAlpha', 0.8);
bg_color = [1, 1, 1];


%% Beautify
axis equal
axis vis3d;
t.EdgeColor = 'none';
set(t,fsoft);
l = light('Position',[-0.2 -0.4 1]);
set(gca,'Visible','off');
set(gcf,'Color',bg_color);
s = add_shadow(t,l,'Color',bg_color*0.8,'BackgroundColor',bg_color,'Fade','infinite');
drawnow


end
