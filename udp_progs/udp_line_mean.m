function udp_line_mean

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
while true
    
    full = fread(u);
    
    if size(full,1) < 8
        continue
    end
    
    iden = arr2int(full(end-3:end),numBytes);
    timestamp = arr2int(full(end-7:end-4), numBytes);
    
    
    if( iden == hex2dec('DEADBEEF') )
        
        for ii=1:4:(size(full,1)-8)
            a = full(ii:ii+3);
            
            rho((ii+3)/4) = arr2int(a,numBytes);
        end
        
        theta = linspace(-3*pi/4,3*pi/4, size(rho,2));
        
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
                
                dx = mean(x - xp);
                dy = mean(y - yp);

                %cla
                hold on
                
                dt = (timestamp-lasttime) * 10^0;
                
                %theta = deg2rad(mod(yaw,360));
                theta = 0;
                R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
                
                
                dx = ([dx dy]*R);
                %dt
                base = base + dx;
                %plot(x-base(1),y-base(2),'b.');
                plot(base(1),base(2), 'go', 'linewidth', 1);
                
                test = [x y]*R;
                x = test(:,1);
                y = test(:,2);
                plot(x-base(1),y-base(2), 'r.');
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

end