
%% distribution of cluster sizes N(n,t)= mean number of clusters of size n at time t
clear N;
N=zeros(700,size(c,2));
T=1:size(c,2);
edges=0.5:1:700.5;
for I=T
            for i=1:size(c,1)
                N(:,I)=N(:,I)+histcounts(c{i,I},edges)';
            end
end
    
I=find(sum(N,2),1,'last');
N=N(1:I,:)/size(c,1);

    %%
    figure(1)
    clf;
    I=2500;
    n=1:size(N,1);
    loglog(n,N(:,I));
    %image(N)



%% mean number of clusters greater than k, number of free particles and size of largest cluster

k=10;
figure(2)
clf;
subplot(3,1,1)
    loglog(t,sum(N(k+1:end,:)))
    hold on;
    loglog(t(1500:2500),1000*(t(1500:2500)).^(-0.25))
    hold off;

subplot(3,1,2)
    semilogx(t,500-[2:size(N,1)]*N(2:end,:))
    ylim([0,500])
    
subplot(3,1,3)
    
    for i=T
    largest(i)=find(N(:,i),1,'last');
    end
    loglog(t,largest)
    
%% mean cluster size bigger than k
k=10;
n=k+1:size(N,1);
tmp=sqrt(n/pi)*N(k+1:end,:);
%tmp=n*N(k+1:end,:);
R=tmp./sum(N(k+1:end,:));

figure(4)
loglog(t,R)
hold on
range=1600:length(t);
power=(4/15);
a=R(range(1))/t(range(1))^(power);
loglog(t(range),a*t(range).^(power))
hold off
%xlim([1,1000])
%ylim([1,10])
xlabel('times (MCS)')
ylabel('linear cluster size')

%%
%alternative method
%mean_size=cellfun(@(x) nanmean(x(x>k)),c);
%figure(5)
%loglog(t,nanmean(mean_size,1))

%xlim([1,100])
%ylim([100,1000])

%%


    