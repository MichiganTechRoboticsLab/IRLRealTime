function [u]=udp_line_icp

u = udp('192.168.0.2', 1140, 'LocalPort', 3540);

numBytes = 4;
set(u,'InputBufferSize',6000)
gg = onCleanup(@() myclean(u));

mag = @(x) sqrt(x(1)^2 + x(2)^2);

fopen(u);

ncorners = 3;

last = [];
skp_ii = 0;
yaw=0;

lasttime = [];

base = [ 0 0 ];

close all
figure
axis([-10 10 -10 10]);

skp = 2;
skpcnt = 0;
while true
    skpcnt = skpcnt + 1;
    full = fread(u);
    
    if size(full,1) < 8 %|| mod(skpcnt,skp) == 0
        continue
    end
    
    iden = arr2int(full(end-3:end),numBytes);
    timestamp = arr2int(full(end-7:end-4), numBytes);
    
    
    if( iden == hex2dec('DEADBEEF') )
        
        for ii=1:4:(size(full,1)-8)
            a = full(ii:ii+3);
            
            rho((ii+3)/4) = arr2int(a,numBytes);
        end
        
        alpha = deg2rad(mod(yaw,360));
        theta = linspace(-3*pi/4,3*pi/4, size(rho,2));
        %theta = mod(theta - alpha, 2*pi);
        
        rho = rho/1000;

        scan_range=rho>1 & rho<15;
        rho = rho(scan_range);
        theta = theta(scan_range);
        %theta = theta(scan_range) - deg2rad(mod(yaw,360));
        
        if isempty(rho)
            continue
        end
        cur = [rho' theta'];
        
        both = [];
        
        if ~isempty(last)
            lasttheta=last(:,2);
            for ii=1:size(cur,1)
                test = find(lasttheta == theta(ii));
                if ~isempty(test)
                    both = [both; cur(ii,:) last(test,:)];
                end
            end
            
            if ~isempty(both)
                x = cos(both(:,2)).*both(:,1);
                y = sin(both(:,2)).*both(:,1);

                xp = cos(both(:,4)).*both(:,3);
                yp = sin(both(:,4)).*both(:,3);
                
                [t, phi] = calcdisp([x y],[xp yp])
                
                dx = real(t(1));
                dy = real(t(2));
                
                if abs(dx) < .003
                    dx = 0;
                else
                    %disp('moving in x')
                    %t(1)
                end
                
                if abs(dy) < .003
                    dy = 0;
                else
                    %disp('moving in y')
                    %t(2)
                end
                if mag([dx dy]) < 0.003 
                    dx=0;dy=0;
                end
                
                %dx = mean(x - xp);
                %dy = mean(y - yp);

                %cla
                hold on
                
                %dt = (timestamp-lasttime) * 10^0;
                
                phi=0;
                R = [cos(phi) -sin(phi);sin(phi) cos(phi)];
                
                
                dxy = ([dx dy]*R);
                dx = dxy(1);
                dy = dxy(2);
                %dt
                base = base + [dx dy];
                dispft = mag(base)*3.28
                
                %plot(x-base(1),y-base(2),'b.');
                
                %cla
                plot(base(1),base(2), 'go', 'linewidth', 3);
                
                %test = [x y]*R;
                %x = test(:,1);
                %y = test(:,2);
                %plot(x-base(1),y-base(2), 'r.');
                axis([-3 3 -3 3]);
                axis square
                %hold off
                pause(0.01);
            end
        end
        
        
        %theta = theta(scan_range)

        last = cur;
        lasttime = timestamp;
        skp_ii = skp_ii+1;
    elseif( iden == hex2dec('CAFEBABE') )
        roll = arr2num(full(1:8),8);
        pitch = arr2num(full(9:16),8);
        yaw = arr2num(full(17:24),8);
    else
        disp('Corrupt')
    end
    
    
    flushinput(u)
    
end

end

function myclean(u)

disp('cleaning up')

fclose(u)
delete(u)

end