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
cn = 10;
stdSet = Inf;
R2min = 0.5;
minldrdist = 0.12;
maxldrdist = .2;
dpmin = 0.3;

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
        
        for ii=1:2:(length(x) - cn)
            
            %%{
            pA = [x(ii) y(ii)];
            pB = [x(ii+floor(cn/2)) y(ii+floor(cn/2))];
            pC = [x(ii+cn-1) y(ii+cn-1)];
            
            pA = pA-pB;
            pC = pC-pB;
            
            dp = dot(pA,pC) / (mag(pA)*mag(pC));
            %}
            
            if( abs(dp) < dpmin && std(pA) < stdSet && std(pC) < stdSet )
                
                [~,~,R2a] = IRLReg(  x(ii:ii+floor(cn/2)),   y(ii:ii+floor(cn/2))  );
                [~,~,R2b] = IRLReg(  x(ii+floor(cn/2):ii+cn-1),   y(ii+floor(cn/2):ii+cn-1)  );
                
                 if( abs(R2a) * abs(R2b) > R2min )   
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
            
            plot(fndpts(:,1),fndpts(:,2),'go', 'linewidth', 3);
            
            keypts = [{fndpts} keypts(1:end-1)];
            if skp_ii > ncorners
                full = [];
                for kk=1:ncorners
                    full = [full; keypts{kk}];
                end
                
                gg = findnclusters(full(:,1:2),1,1,minldrdist);
                
                plot(gg(:,1),gg(:,2),'mo', 'linewidth', 4);
                
                calcmean = [mean(gg(:,1)) mean(gg(:,2))];
                
                %plot(calcmean(1), calcmean(2), 'kx', 'linewidth', 5);
            end
            %{
            if(~isempty(ctrs) && size(ctrs,2) > 1)
                corners = [{ctrs} corners(1:end-1)]
                
                if( skp_ii > ncorners )
                    
                    maxn = 0;
                    full = [];
                    
                    for jj=1:ncorners-1
                        
                        full = [full; corners{jj}];
                        corn_n = size(corners{jj},1);
                        if maxn < corn_n
                            maxn = corn_n;
                        end
                        
                        temp = corners{jj};
                        plot(temp(:,1),temp(:,2),'r+', 'linewidth', 3);
                    end
                    
                    %{
                    for kk=1:size(full,1)
                        tmpfull = full;
                        kk
                        tmpfull(kk,:) = [];
                        
                        full(kk,:)
                        ret = findclosept(full(kk,:), tmpfull, minldrdist);
                        
                        if ret == -1
                            full = tmpfull;
                        end
                        
                        if kk>= size(full)
                            break
                        end
                    end
                    %}
                    [idx,ctrs_adj] = kmeans(full(:,1:2), maxn, 'emptyaction', 'drop');
                    
                    plot(ctrs_adj(:,1),ctrs_adj(:,2),'go', 'linewidth', 2);
                    
                end
                
                %{
                if( skp_ii >= ncorners && ~isempty(corners{ncorners}) )
                    maxn = 0;
                    full = [];
                    
                    for jj=1:ncorners
                        
                        full = [full; corners{jj}];
                        corn_n = size(corners{jj},1);
                        if maxn < corn_n
                            maxn = corn_n;
                        end
                        
                        temp = corners{jj};
                        plot(temp(:,1),temp(:,2),'r+', 'linewidth', 3);
                    end
                    
                    maxn
                    full
                    [idx,~] = kmeans(full(:,1:2), maxn, 'emptyaction', 'drop');
                    
                    
                    tokeep = zeros(length(idx),1);
                    shiftx = 0;
                    shifty = 0;
                    
                    for kk=1:max(idx)
                        
                        tokeep(kk) = length(idx(idx==kk)) > floor(ncorners*.7);
                        
                        if tokeep(kk);
                            
                            a = find(idx == kk);
                            
                            tempx = full(a(end),1)-full(a(1),1);
                            shiftx = [shiftx tempx];
                            
                            tempy = full(a(end),2)-full(a(1),2);
                            shifty = [shifty tempy];
                            
                        end
                        
                        
                        
                    end
                    
                    shiftx = shiftx(abs(shiftx) < 1.5*std(shiftx));
                    shifty = shifty(abs(shifty) < 1.5*std(shifty));
                    
                    shift = [mean(shiftx) mean(shifty)];
                    
                    
                    if shift(1) < minldrdist || shift(1) > maxldrdist
                        shift(1) = 0;
                    end
                    
                    if shift(2) < minldrdist || shift(2) > maxldrdist
                        shift(2) = 0;
                    end
                    curcenter = shift + curcenter
                    
                    for jj=1:length(corners)
                        
                        temp = corners{jj};
                        for kk=1:size(temp,1)
                            
                            temp(kk,:) = temp(kk,:) - shift;
                            
                        end
                        corners{jj} = temp;
                        
                    end
                end
                %}
            end
            %}
            
            skp_ii = skp_ii+1;
        end
        
        hold off
        
        
        pause(.01)
    elseif( iden == hex2dec('CAFEBABE') )
        roll = arr2num(full(1:8),8);
        pitch = arr2num(full(9:16),8);
        yaw = arr2num(full(17:24),8);
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