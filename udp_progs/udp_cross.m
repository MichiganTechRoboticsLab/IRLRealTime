function udp_cross

u = udp('192.168.0.2', 1140, 'LocalPort', 3540);
numBytes = 4;
set(u,'InputBufferSize',6000)
ret = onCleanup(@() myclean(u));

skp_ind= 1;
skp_cnt= 10;

fopen(u);
yaw=0;

while true
        
    full = fread(u);
    
    if size(full,1) < 8
        continue
    end
    
    iden = arr2int(full(end-3:end),numBytes);
    timestamp = arr2int(full(end-7:end-4), numBytes);
    
    
    hold on
    if( iden == hex2dec('DEADBEEF') )
        for ii=1:4:(size(full,1)-8)
            a = full(ii:ii+3);
            
            rho((ii+3)/4) = arr2int(a,numBytes);
        end
        
        %size(full)
        
        theta = linspace(-3*pi/4,3*pi/4, size(rho,2));
        
        rho = rho/1000;
        
        scan_range=rho>1.5 & rho<10;
        rho = rho(scan_range);
        length(rho)
        
        if( isempty(rho) )
            continue
        end
        
        theta = theta(scan_range) - deg2rad(mod(yaw,360));
        
        x = cos(theta).*rho;
        y = sin(theta).*rho;
        
        cla
        %axis([-4 8 -2 3])
        axis([-4 4 -4 4])
        %axis([-0.5 0.5 -0.5 0.5])
        
        %if( mod(skp_ind,skp_cnt) == 0 )
            plot(x,y, 'b.')
        %end
        
        skp_ind = skp_ind + 1;
        
        hold off
        pause(.001)
    elseif( iden == hex2dec('CAFEBABE') )
        %disp('found vectornav')
        
        roll = arr2num(full(1:8),8);
        pitch = arr2num(full(9:16),8);
        yaw = arr2num(full(17:24),8);
        
        %[roll pitch mod(yaw,360)]
    else
        %disp('Corrupt')
        %size(full)
    end
    
    flushinput(u)
    
end
fclose(u)

end

function myclean(u)

disp('cleaning up')

fclose(u)

end