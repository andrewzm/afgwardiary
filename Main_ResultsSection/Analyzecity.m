function Analyzecity(s1star,s2star,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)

% Plot insets on growth map figure

% Initialize
i=1;
k=1;
Weeknum = zeros(floor(length(spikeAll)/7),1);
N = length(Weeknum);
subplot(1,2,1)
hold on

% Organise data once again into weeks (this is redundant...) and scatter
% events on map
for k = 1:floor(length(spikeAll)/7)
    for j = 0:6  %Cluster in weeks
        if ~isempty(spikeAll(i).Coords)
            for l = 1:size((spikeAll(i).Coords),1)
                if (spikeAll(i).Coords(l,1) < ucorner(1)) && (spikeAll(i).Coords(l,1) > lcorner(1)) && (spikeAll(i).Coords(l,2) < ucorner(2)) && (spikeAll(i).Coords(l,2) > lcorner(2))
                    Weeknum(k) = Weeknum(k) + 1;
                end
            end
            plot(spikeAll(i).Coords(:,1),spikeAll(i).Coords(:,2),'k.')
        end
        i = i+1;
        hold on
    end
end
xlabel('Lon','FontSize',14)
ylabel('Lat','FontSize',14)
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'FontSize',14)
axis([lcorner(1) ucorner(1) lcorner(2) ucorner(2) ])
box off

% Now we can plot the smoothed intensity
subplot(1,2,2)
plot(Weeknum,'k')

lon1 = find(abs(s1star - lcorner(1)) == min(abs(s1star - lcorner(1))));
lon2 = find(abs(s1star - ucorner(1)) == min(abs(s1star - ucorner(1))));
lat1 = find(abs(s2star - lcorner(2)) == min(abs(s2star - lcorner(2))));
lat2 = find(abs(s2star - ucorner(2)) == min(abs(s2star - ucorner(2))));

hold on
medianintensityest = zeros(N,1);
upperintensityest = zeros(N,1);
lowerintensityest = zeros(N,1);

for i = 1:N
    medianintensityest(i) = sum(sum(int_median(lat1:lat2,lon1:lon2,i)))*ds^2; 
    upperintensityest(i) = sum(sum(int_upper(lat1:lat2,lon1:lon2,i)))*ds^2; 
    lowerintensityest(i) = sum(sum(int_lower(lat1:lat2,lon1:lon2,i)))*ds^2; 
end

jbfill([1:N-1],upperintensityest(1:N-1)',lowerintensityest(1:N-1)','g','g',1,0.5)
hold on;

% Labels
xlabel('k','FontSize',14)
ylabel('N_S','FontSize',14)
axis('tight')
set(gca,'Xlim',[0 320])
set(gca,'XTick',[0 320])
set(gca,'FontSize',14)
box off

