SC = importdata("SC.dat");
Mars = importdata("Mars.dat");
Earth = importdata("Earth.dat");

figure
plot(Earth(:,2), Earth(:,3))%, Earth(:,4));
hold on
plot(Mars(:,2), Mars(:,3))%, Mars(:,4));
plot(SC(:,2), SC(:,3))%, SC(:,4));
axis equal