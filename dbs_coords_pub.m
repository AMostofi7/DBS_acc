clear all

[database path] = uigetfile('*.*'); % Defines the database file (see database template excel file) and path

sheet = 'Sheet3'; % Relevant sheet of the database

[dets, heads] = xlsread([path, database], sheet); % Reads excel database. "Dets" matrix contains numerical values from spreadsheet.

% TRAJECTORY VECTORS AND CONTACT COORDINATES

dets2 = dets(~isnan(dets(:,3)),:);
dets2 = dets2(:,[1, 3:end]);

planned = dets(:,[19:25]); % [target x,y,z entry x,y,z]
actual = dets(:,[26:32]); % [tip x,y,z shaft x,y,z]

for n=1:size(planned,1)
plan_unit(n,:) = [planned(n,1:3) - planned(n,4:6)] ./ norm([planned(n,1:3) - planned(n,4:6)]);
act_unit(n,:) = [actual(n,1:3) - actual(n,4:6)] ./ norm([actual(n,1:3) - actual(n,4:6)]);
end

planned_tip = planned(:,1:3) + planned(:,7).*plan_unit; % add unit*offset to target

error3D = [actual(:,1:3) - planned_tip, dets(:,1)]; 
error3DR = error3D(find(planned(:,1)>0),:);
error3DL = error3D(find(planned(:,1)<0),:);

% 3D scatter plot of actual electrode tips from intended tip (origin)
figure
scatter3(error3DR(:,1), error3DR(:,2), error3DR(:,3),'b');
hold on
scatter3(error3DL(:,1), error3DL(:,2), error3DL(:,3),'r');
xlabel('X')
ylabel('Y')
zlabel('Z')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.ZAxisLocation = 'origin';

tgt_refs = (actual(:,3) - planned(:,3))./act_unit(:,3); % distance between actual tip Z and planned target Z in unit Z components
short_refs = find(tgt_refs < 0); % electrodes short of the target
                                % tip minus unit Z-component
tgt_refs(short_refs,:) = NaN;                                
act_td = actual(:,1:3) - act_unit.*tgt_refs; % actual at target depth
est_td_shorts = NaN(size(act_td,1), size(act_td,2)); % estimate for those short electrodes
est_td_shorts(short_refs,:) = actual(short_refs,1:3); % use tips as best estimate for short electrodes

errorXY = act_td(:,1:2) - planned(:,1:2);
errorXYshorts = est_td_shorts(:,1:2) - planned(:,1:2);

errorXYR = errorXY(find(planned(:,1)>0),:);
errorXYL = errorXY(find(planned(:,1)<0),:);

figure
scatter(errorXYR(:,1), errorXYR(:,2), 'bo');
%scatter(errorXYR(end-20:end,1), errorXYR(end-20:end,2), 'bo');
hold on
scatter(errorXYL(:,1), errorXYL(:,2), 'ro');
%scatter(errorXYL(end-20:end,1), errorXYL(end-20:end,2), 'ro');

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
viscircles([ 0 0; 0 0; 0 0], [1,2,3], 'color' ,'k', 'LineStyle', ':', 'Linewidth', 0.5 )
grid on
xlabel('X')
ylabel('Y')

% ERROR BARS (point X, point Y, yneg, ypos, xneg, xpos, marker)
errorbar(nanmean(errorXYR(:,1)), nanmean(errorXYR(:,2)), nanstd(errorXYR(:,2)), nanstd(errorXYR(:,2)), nanstd(errorXYR(:,1)), nanstd(errorXYR(:,1)), 'db')
errorbar(nanmean(errorXYL(:,1)), nanmean(errorXYL(:,2)), nanstd(errorXYL(:,2)), nanstd(errorXYL(:,2)), nanstd(errorXYL(:,1)), nanstd(errorXYL(:,1)), 'dr')


% Trajectories [hosp no, target X, target Y, target Z, vector X, vector Y, vector Z] 
rt_plan = [dets2(:,[1, 2:4]), normr([(dets2(:,2:4)-dets2(:,5:7))])];
lt_plan = [dets2(:,[1, 8:10]), normr([(dets2(:,8:10)-dets2(:,11:13))])];
rt_act = [dets2(:,[1, 14:16]), normr([(dets2(:,14:16)-dets2(:,17:19))])];
lt_act = [dets2(:,[1, 20:22]), normr([(dets2(:,20:22)-dets2(:,23:25))])];

% Contacts [X0, Y0, Z0, X1, Y1, Z1, ...]
rt_plan_conts = [rt_plan(:,2:4)-cont(1,2)*rt_plan(:,5:7), rt_plan(:,2:4)-cont(2,2)*rt_plan(:,5:7), rt_plan(:,2:4)-cont(3,2)*rt_plan(:,5:7), rt_plan(:,2:4)-cont(4,2)*rt_plan(:,5:7)];
rt_act_conts = [rt_act(:,2:4)-cont(1,2)*rt_act(:,5:7), rt_act(:,2:4)-cont(2,2)*rt_act(:,5:7), rt_act(:,2:4)-cont(3,2)*rt_act(:,5:7), rt_act(:,2:4)-cont(4,2)*rt_act(:,5:7)];

lt_plan_conts = [lt_plan(:,2:4)-cont(1,2)*lt_plan(:,5:7), lt_plan(:,2:4)-cont(2,2)*lt_plan(:,5:7), lt_plan(:,2:4)-cont(3,2)*lt_plan(:,5:7), lt_plan(:,2:4)-cont(4,2)*lt_plan(:,5:7)];
lt_act_conts = [lt_act(:,2:4)-cont(1,2)*lt_act(:,5:7), lt_act(:,2:4)-cont(2,2)*lt_act(:,5:7), lt_act(:,2:4)-cont(3,2)*lt_act(:,5:7), lt_act(:,2:4)-cont(4,2)*lt_act(:,5:7)];

% Error
for n = 1:size(dets2,1)
    
    %SCALAR 
%     % error reference, midpoint of contacts between 1 and 2
%     rt_erref = 0.5*(rt_plan_conts(n, 4:6) + rt_plan_conts(n, 7:9));
%     lt_erref = 0.5*(lt_plan_conts(n, 4:6) + lt_plan_conts(n, 7:9));
    
    % actual X,Y difference from target Z / unit Z vector   X,Y unit vectors 
   rt_ref = [ (rt_act(n,2:3) + ((rt_plan(n,4) - rt_act(n,4))/rt_act(n,7))*rt_act(n,5:6)) rt_plan(n,4)];
   
   rt_act_td = rt_act(n,2:3) + ((rt_plan(n,4) - rt_act(n,4))/rt_act(n,7))*rt_act(n,5:6); % X-Y coords at plan target Z
   lt_act_td = lt_act(n,2:3) + ((lt_plan(n,4) - lt_act(n,4))/lt_act(n,7))*lt_act(n,5:6);
   
   rt_scalar = rt_act_td - rt_plan(n,2:3);
   lt_scalar =  lt_act_td - lt_plan(n,2:3);
    
     % RADIAL - arccos(dot product)
    
    rt_error(n,:) = [dets(n,1), rt_scalar, norm(rt_scalar),  acosd(dot(rt_act(n,5:7), rt_plan(n,5:7)))];
    lt_error(n,:) = [dets(n,1), lt_scalar, norm(lt_scalar),  acosd(dot(lt_act(n,5:7), lt_plan(n,5:7)))];
    
end

figure
scatter(rt_error(:,2), rt_error(:,3),'bo')
hold on
scatter(lt_error(:,2), lt_error(:,3),'ro')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

plot(mean(rt_error(:,2)), mean(rt_error(:,3)), 'db', 'MarkerFaceColor', 'b', 'Markersize', 10)
plot(mean(lt_error(:,2)), mean(lt_error(:,3)), 'dr', 'MarkerFaceColor', 'r', 'Markersize', 10)

viscircles([ 0 0; 0 0], [1.27 , 2.54], 'color' ,'k', 'Linewidth', 0.5 )

% ACTIVE CONTACT COORDINATES
rt_AC_coords = dets2(:,[1, 26, 28:30]);
lt_AC_coords = dets2(:,[1, 27, 31:33]);

lt_stim = dets2(:,[1 41:46]);
rt_stim = dets2(:,[1 35:40]);

% ACTIVE CONTACT COORDINATES FOR EACH ANATOMICAL GROUP
[gdets, gheads] = xlsread([path database],'Sheet6');
gdets = gdets(:,1:8);

clear ac_grp stim_grp
r = 2 %change for different groups (1-4)
for n=1:size(gdets,1)
    
    if ~isnan(gdets(n,(2*r-1)))
    if gdets(n,2*r) == 0
        ac_grp(n,:) = rt_AC_coords(find(rt_AC_coords(:,1) == gdets(n,(2*r-1))), 2:5);
        ac_grp(n,1) = ac_grp(n,1)-8;
        
        stim_grp(n,:) = rt_stim(find(rt_stim(:,1) == gdets(n,(2*r-1))), 2:end);
        
    elseif gdets(n,2*r) == 1
        ac_grp(n,:) = lt_AC_coords(find(lt_AC_coords(:,1) == gdets(n,(2*r-1))), 2:5);
    stim_grp(n,:) = lt_stim(find(lt_stim(:,1) == gdets(n,(2*r-1))), 2:end);
    end
    end
end

ac_grp2 = [abs(ac_grp(:,2)) ac_grp(:,3:4)]; % MODULUS OF X-coordinate
% xlswrite([path database],ac_grp2,'Sheet14','j1');

mean(ac_grp2)
std(ac_grp2)/sqrt(size(ac_grp,1));
