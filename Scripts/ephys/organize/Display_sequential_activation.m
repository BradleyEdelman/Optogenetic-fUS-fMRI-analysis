load('D:\Data_Processed\ephys\20191217_4364123_R1\1_0\1_0.mat')
%%
figure; plot(LFP_mV{1})

%%
t_bank = {'L CPu', 'L Cortex', 'R CPu', 'R Cortex'};

C = [253,184,99; 230,97,1; 178,171,210; 94,60,153]/255;

t1 = 3.25e5;
t2 = 5e5;

l1 = LFP_mV{1}(t1:t2)*2-2.5;
l2 = LFP_mV{16}(t1:t2)*2-2;
l3 = LFP_mV{17}(t1:t2)-1;
l4 = LFP_mV{32}(t1:t2);

figure(1); clf; hold on
plot(1:size(l1,1),l1,'color',C(1,:));
plot(1:size(l1,1),l2,'color',C(2,:));
plot(1:size(l1,1),l3,'color',C(3,:));
plot(1:size(l1,1),l4,'color',C(4,:));
set(gca,'ylim',[-3.25 .5])

figure(2); clf; hold on
plot(1:size(l1,1),l1,'color',C(1,:));
plot(1:size(l1,1),l2,'color',C(2,:));
plot(1:size(l1,1),l3,'color',C(3,:));
plot(1:size(l1,1),l4,'color',C(4,:));
set(gca,'ylim',[-3.25 .5],'xlim',[8.5e4 10e4])

figure(3); clf; hold on
plot(1:size(l1,1),l1,'color',C(1,:),'linewidth',2);
plot(1:size(l1,1),l2,'color',C(2,:),'linewidth',2);
plot(1:size(l1,1),l3,'color',C(3,:),'linewidth',2);
plot(1:size(l1,1),l4,'color',C(4,:),'linewidth',2);
t11 = 9.035e4;
t21 = 9.065e4;
set(gca,'ylim',[-3.25 .5],'xlim',[t11 t21])


t1 = [t11:1e2:t21];
t2 = [t21-2e2:1e2:t21+1e2];
for i = 1:size(t1,2)
    tt = t1(i):t2(i);
    l11 = l1(tt); xline(tt(find(l11 == min(l11))),'color',C(1,:),'linewidth',1.5)
    l21 = l2(tt); xline(tt(find(l21 == min(l21))),'color',C(2,:),'linewidth',1.5)
    l31 = l3(tt); xline(tt(find(l31 == min(l31))),'color',C(3,:),'linewidth',1.5)
    l41 = l4(tt); xline(tt(find(l41 == min(l41))),'color',C(4,:),'linewidth',1.5)
end

% set(gca,'xlim',[9.045 9.055]*10e3)


%%
t1 = [5e4 1.1e5 1.7e5 2.3e5 2.9e5 3.e5 4.1e5 4.7e5 5.3e5 5.9e5];
t2 = [7e4 1.25e5 1.84e5 2.45e5 3.05e5 3.65e5 4.25e5 4.85e5 5.45e5 6.05e5];
for i = 1:10
    clear dtmp
    dtmp(1,:) = LFP_mV{1}(t1(i):t2(i));
    dtmp(2,:) = LFP_mV{16}(t1(i):t2(i));
    dtmp(3,:) = LFP_mV{17}(t1(i):t2(i));
    dtmp(4,:) = LFP_mV{32}(t1(i):t2(i));
    
    [pks,locs] = findpeaks(var(dtmp,1));
    pks2 = find(pks>0.1); pks2 = pks2(1); %first real peak
    
    t11 = locs(pks2) - 10;
    t21 = locs(pks2) + 10;
    
    for ii = 1:120-1
        t11 = t11 + 100;
        t21 = t21 + 100;
        dtmp2 = dtmp(:,t11:t21);
        for j = 1:4
            tmp = dtmp2(j,:);
            acttmp = find(tmp == min(tmp));
            act(j,ii,i) = acttmp(1);
        end
    end
    
    act(:,:,i) = act(:,:,i) - repmat(act(end,:,i),4,1);
    
end
act = reshape(act,4,[]);
%%

figure; f = bar(flipud(mean(act,2))); set(gca,'xticklabel',fliplr(t_bank))
hold on
errorbar(flipud(mean(act,2)),flipud(std(act,[],2)))
    
