function finalplot(param, phases, res)

loadconstants


colormatrix = [ 0         0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
colindex = 1;
for ip = 1:numel(phases)
    colour_found = 0;
    if ip > 1
        for ifind = ip-1:-1:1
            if isequal(phases(ip).vehicle,phases(ifind).vehicle)
                phases(ip).plotcolor = phases(ifind).plotcolor;
                colour_found = 1;
            end
        end
    end
    if colour_found == 0;
        phases(ip).plotcolor = colormatrix(colindex,:);
        colindex = colindex + 1;
    end
end

% add bounds to figure 1
figure(1)
index = 0;
for ip = 1:numel(phases)
    for ine = 1:phases(ip).ne
        index=index+1;
        subplot(4,2,1),hold on
        xlabel('Time [s]')
        ylabel('AoA [deg]')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).cbounds(1,1), phases(ip).cbounds(1,1)]),'k:')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).cbounds(1,2), phases(ip).cbounds(1,2)]),'k:')
        
        subplot(4,2,2),hold on
        xlabel('Time [s]')
        ylabel('Throttle')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).cbounds(2,1), phases(ip).cbounds(2,1)],'k:')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).cbounds(2,2), phases(ip).cbounds(2,2)],'k:')
        
        subplot(4,2,3),hold on
        xlabel('Time [s]')
        ylabel('Roll [deg]')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).cbounds(3,1), phases(ip).cbounds(3,1)]),'k:')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).cbounds(3,2), phases(ip).cbounds(3,2)]),'k:')

        subplot(4,2,4),hold on
        xlabel('Time [s]')
        ylabel('Altitude [km]')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(1,1)/1e3, phases(ip).xbounds(1,1)/1e3],'k:')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(1,2)/1e3, phases(ip).xbounds(1,2)/1e3],'k:')

        subplot(4,2,5),hold on
        xlabel('Time [s]')
        ylabel('Relative velocity [km/s]')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(2,1)/1e3, phases(ip).xbounds(2,1)/1e3],'k:')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(2,2)/1e3, phases(ip).xbounds(2,2)/1e3],'k:')

        subplot(4,2,6),hold on
        xlabel('Time [s]')
        ylabel('Relative FPA [deg]')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).xbounds(3,1), phases(ip).xbounds(3,1)]),'k:')
        plot([res(index).t(1), res(index).t(end)], rad2deg([phases(ip).xbounds(3,2), phases(ip).xbounds(3,2)]),'k:')
        
        subplot(4,2,7),hold on
        if ip==1 && ine==1
            geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
            minlat = +pi;   maxlat = -pi;   minlon = pi/2;  maxlon = -pi/2;
        end
        minlat = min(min(res(index).x(:,5)),minlat);
        maxlat = max(max(res(index).x(:,5)),maxlat);
        minlon = min(min(res(index).x(:,6)),minlon);
        maxlon = max(max(res(index).x(:,6)),maxlon);
        if ip==numel(phases) && ine==phases(ip).ne
            axis(rad2deg([minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1]))
        end
        
        subplot(4,2,8),hold on
        xlabel('Time [s]')
        ylabel('Mass [ton]')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(7,1)/1e3, phases(ip).xbounds(7,1)/1e3],'k:')
        plot([res(index).t(1), res(index).t(end)], [phases(ip).xbounds(7,2)/1e3, phases(ip).xbounds(7,2)/1e3],'k:')
    end
end

%% add bounds to figure 2
figure (2)
subplot(2,3,4),hold on
plot([0 res(end).t(end)],[1 1], 'k:')


%% create other figures
E = wgs84Ellipsoid;

index=0;
for ip = 1:numel(phases)
    for ine = 1:phases(ip).ne
        index=index+1;
        
        figure(3)
        [xeci,yeci,zeci]=enu2ecef(0,0,0,res(index).x(:,5),res(index).x(:,6),res(index).x(:,1),E,'radians');
        plot3(xeci,-yeci,zeci, 'color', phases(ip).plotcolor, 'LineWidth', 2), hold on
        if isfield(param,'no_fly')
            if ~isempty(param.no_fly)
                no_fly_arc = param.no_fly(:,3)./const.rE;
                [nfzlat,nfzlon] = scircle1(param.no_fly(:,1),param.no_fly(:,2),no_fly_arc);
                theta=-pi:pi/100:pi;
                xcirc = cos(theta);
                ycirc = sin(theta);
                zcirc = 0*theta;
                for ici = 1 : size(nfzlat,2)
                    [xnfz,ynfz,znfz] =  enu2ecef(param.no_fly(ici,3).*xcirc,param.no_fly(ici,3).*ycirc,zcirc,...
                                        param.no_fly(ici,1),param.no_fly(ici,2),const.rE*cos(no_fly_arc(ici)),...
                                        referenceSphere,'radians');
                    plot3(xnfz,-ynfz,znfz, 'color', 'r', 'LineWidth', 1)
                end
            end
        end
        
        figure(4)
        plot3(wrapTo180(rad2deg(res(index).x(:,5))), wrapTo180(rad2deg(res(index).x(:,6))),...
            res(index).x(:,1).*1e-3,'Color',phases(ip).plotcolor), hold on
        
        figure(5)
            segment_distances = distance_haversine (res(1).x(1,5), res(1).x(1,6),...
                res(index).x(:,5), res(index).x(:,6), const.rE);
        plot(segment_distances./1e3,res(index).x(:,1)./1e3,'Color',phases(ip).plotcolor), hold on
    end
end
%%
figure(3)
[xnfz,ynfz,znfz]=sphere(100);
Xt=const.rE*xnfz;
Yt=const.rE*ynfz;
Zt=const.rE*znfz;
rotate3d
I=imread('earth_blb.png');
warp(Xt,Yt,Zt,I),        hold on
axis equal
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
view(90+(rad2deg(res(1).x(1,6))+rad2deg(res(end).x(end,6)))/2,...
    (rad2deg(res(1).x(1,5))+rad2deg(res(end).x(end,5)))/2)
%%
figure(4)
view(45,45)
axis tight
xm = mean(xlim); ym = mean(ylim); % find average value
halfx = abs(diff(xlim))/2;  halfy = abs(diff(ylim))/2; % get size
hs=max(halfx,halfy); % find max
xlim ([xm-hs,xm+hs]); ylim ([ym-hs,ym+hs]); %apply new limits
for index = 1:numel(res)
    plot3((xm+hs).*ones(size(res(index).x(:,5))), wrapTo180(rad2deg(res(index).x(:,6))),...
                (res(index).x(:,1)).*1e-3, 'Color', 0.6.*[1 1 1])
    plot3(wrapTo180(rad2deg(res(index).x(:,5))), (ym+hs).*ones(size(res(index).x(:,6))),...
                (res(index).x(:,1)).*1e-3, 'Color', 0.6.*[1 1 1]) 
end
xlabel('Latitude (deg)'), ylabel('Longitude (deg)'), zlabel('Altitude (km)')
box on     
%%
figure(5)
xlabel('Downrange Distance (km)'), ylabel('altitude (km)')       