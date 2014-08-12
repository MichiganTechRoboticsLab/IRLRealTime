function udp_line

u = udp('192.168.0.2', 1140, 'LocalPort', 3540);
numBytes = 4;
set(u,'InputBufferSize',6000)
ret = onCleanup(@() myclean(u));

mag = @(x) sqrt(x(1)^2 + x(2)^2);

fopen(u);

ncorners = 3;
corners = cell(1,ncorners);
skp_ii = 0;
yaw=0;
cn = 15;
stdSet = 0.3;
R2min = 0.5;
minldrdist = .3;
maxldrdist = .2;
dpmin = 0.6;

keypts = cell(1,3);
curcenter = [0 0];

lastmean = [];

calc_mean_still_n = 10;
mean_still = [];
mean_still_vec = [];

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

        scan_range=rho>1 & rho<10;
        rho = rho(scan_range);
        
        if isempty(rho)
            continue
        end
        theta = theta(scan_range) - deg2rad(mod(yaw,360));
        %theta = theta(scan_range)
        x = cos(theta).*rho;
        y = sin(theta).*rho;
        
        cla
        %axis([-0.5 4 -3.2 0.75])
        %axis([-6 6 -6 6])
        %axis([-0.5 0.5 -0.5 0.5])
        axis square
        hold on
        plot(x,y, 'b.')
        %plot(curcenter(1),curcenter(2),'go','linewidth', 5);
        
        %try some corner detection
        temp = [];
        
        estcluster = 1;
        last = 0;
        
        for ii=1:1:(length(x) - cn)
            
            %%{
            indA = ii;
            indB = ii+floor(cn/2);
            indC = ii+cn-1;
            
            pA = [x(indA) y(indA)];
            pB = [x(indB) y(indB)];
            pC = [x(indC) y(indC)];
            
            pA = pA-pB;
            pC = pC-pB;
            
            dp = dot(pA,pC) / (mag(pA)*mag(pC));
            %}
            
            if( abs(dp) < dpmin && std(pA) < stdSet && std(pC) < stdSet )
                
                la = [x(indA:indB),y(indA:indB)];
                lb = [x(indB:indC),y(indB:indC)];
                
                [~,~,R2a] = IRLReg(x(indA:indB),y(indA:indB));
                [~,~,R2b] = IRLReg(x(indB:indC),y(indB:indC));
                
                 if( abs(R2a) * abs(R2b) * cov(la, lb) > R2min) 
                    if ~last
                        estcluster = estcluster + 1;
                    end
                    
                    %plot(x(ii:ii+cn-1), y(ii:ii+cn-1), 'r.');
                    %plot(pB(1), pB(2), 'mx', 'linewidth', 6);
                    temp = [temp; [pB(1) pB(2)]];
                    last = 1;
                    
                else
                    last = 0;
                end
            else
                last = 0;
            end
        end
        
        
        if estcluster > size(temp,1)
            estcluster = size(temp,1);
        end
        
        if ~isempty( temp )
            gg = findnclusters(temp(:,1:2),1,1,minldrdist);
            estcluster = max(gg(:,3));
            
            %[idx,ctrs,sumd] = kmeans(temp(:,1:2), estcluster, 'emptyaction', 'drop');
            
            plot(temp(:,1),temp(:,2),'ro', 'linewidth', 2);
            
            fndpts = [];
            for jj=1:estcluster
                
                alpha = gg(gg(:,3)==jj,1:2);
                fndpts = [fndpts;mean(alpha(:,1)) mean(alpha(:,2))];
                
            end
            
            %plot(fndpts(:,1),fndpts(:,2),'go', 'linewidth', 3);
            
            keypts = [{fndpts} keypts(1:end-1)];
            if skp_ii > ncorners
                full = [];
                for kk=1:ncorners
                    full = [full; keypts{kk}];
                end
                
                gg = findnclusters(full(:,1:2),1,1,minldrdist);
                
                
                bb = [];
                size(gg)
                for ii=1:max(gg(:,3))
                    temp = gg(gg(:,3) == ii,:);
                    bb = [bb; mean(temp(:,1)) mean(temp(:,2))];
                end
                plot(bb(:,1),bb(:,2),'mo', 'linewidth', 4);
                size(bb)
                %calcmean = [mean(gg(:,1)) mean(gg(:,2))];
                
                %plot(calcmean(1), calcmean(2), 'kx', 'linewidth', 5);
            end
            
            skp_ii = skp_ii+1;
        end
        
        hold off
        
        
        pause(.01)
    elseif( iden == hex2dec('CAFEBABE') )
        roll = arr2num(full(1:8),8) ;
        pitch = arr2num(full(9:16),8);
        yaw = arr2num(full(17:24),8) + 90;
    else
        disp('Corrupt')
    end
    
    
    flushinput(u)
    
end
fclose(u)

end

function myclean(u)

disp('cleaning up')

fclose(u)

end