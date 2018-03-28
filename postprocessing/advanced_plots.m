function advanced_plots (phases, res, forces, accelerations, Mn, thermal)

loadconstants

% initialize plot style
plotcol = 'k';
% lsAB = ':';        lsR = '-';
plotmar = 'none';
plotlw  = 1;
plotst = '-';
% 
% plotT = 1;
% plotACC = 1;
% plotMn =1;

if ishandle(2) && strcmp(get(2, 'type'), 'figure')  % if figure 2 exists
    a = get(2);
    b = get(a.Children(:));
    for ifi = 1:size(b,1)                  % dim the color of the subplots
        set(b(ifi).Children(:),'Color',[0.8 0.8 0.8])
    end
    figure (2)
else                                                % else create it
    scrsz = get(groot,'ScreenSize');
    figure('Color',[1 1 1],'Position',[scrsz(3)/2 50 scrsz(3)/2-10 scrsz(4)-135]);
end

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


index = 0;
con_counter = 0;
for ip = 1:numel(phases)
    for ine = 1:phases(ip).ne
        index=index+1;

        time_points = size(res(index).t,2);
        
        L = forces(con_counter+1:con_counter+time_points,1);
        D = forces(con_counter+1:con_counter+time_points,2);
        FT = forces(con_counter+1:con_counter+time_points,3);
        q = forces(con_counter+1:con_counter+time_points,4);        
        XZacc = accelerations(con_counter+1:con_counter+time_points,:);        
        acc = [XZacc, hypot(XZacc(:,1),XZacc(:,2))];
        Mach = Mn(con_counter+1:con_counter+time_points);   
        temp = thermal(con_counter+1:con_counter+time_points);   

        con_counter = con_counter + time_points;
        
        subplot(2,3,1),hold on
        xlabel('Time [s]')
        ylabel('Aerodyn. Efficiency')
        plot (res(index).t,L./D,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),L(end)/D(end),'k+')
        
        subplot(2,3,2),hold on
        xlabel('Time [s]')
        ylabel('Mach number')
        plot (res(index).t,Mach,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),Mach(end),'k+')
        
        subplot(2,3,3),hold on
        xlabel('Time [s]')
        ylabel('forces [kN]')
        plot (res(index).t,FT./1e3,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),FT(end)./1e3,'k+')
        plot (res(index).t,L./1e3,...
            'Color', 'k','Marker',plotmar,'LineWidth',plotlw-0.5,'LineStyle',plotst)
        plot(res(index).t(end),L(end)./1e3,'k+')        
        plot (res(index).t,D./1e3,...
            'Color', 'k','Marker',plotmar,'LineWidth',plotlw-0.5,'LineStyle',plotst)
        plot(res(index).t(end),D(end)./1e3,'k+')
        
        subplot(2,3,4),hold on
        xlabel('Time [s]')
        ylabel('Acceleration/g0')
        plot (res(index).t,acc(:,1)./9.81,res(index).t,acc(:,2)./9.81,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot (res(index).t,acc(:,3)./9.81,...
            'Color', 'k','Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),acc(end,1)./9.81,'k+')
        plot(res(index).t(end),acc(end,2)./9.81,'k+')  
        plot(res(index).t(end),acc(end,3)./9.81,'k+')
        
        subplot(2,3,5),hold on
        xlabel('Time [s]')
        ylabel('Stagnation temp. [K]')
        plot (res(index).t,temp,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),temp(end,1),'k+')
        
        subplot(2,3,6),hold on
        xlabel('Time [s]')
        ylabel('Dynamic pressure [Pa]')
        plot (res(index).t,q,...
            'Color', phases(ip).plotcolor,'Marker',plotmar,'LineWidth',plotlw,'LineStyle',plotst)
        plot(res(index).t(end),q(end),'k+')
    end
end