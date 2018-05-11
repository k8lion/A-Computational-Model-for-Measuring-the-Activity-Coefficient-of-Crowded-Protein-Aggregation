% load('numsamples2.mat');
% load('numbonds2.mat');
% load('onrateconstants2.mat');
% curr = size(y,1);
load('randwalk_output_test12.0.mat');
newttB = [];
newc = [];
newtime = [];
trials = 12;

parfor i = 1:trials
    ts = tic;
    [times,locs,types,rs,orients,vf1,tt,tD,eD,bP,nB,newttB(i,:),newc(i),eR,aR] = RandWalkLM2(18,18,5,2,0,1000,1000,0,5);
    newtime(i) = toc(ts);
end

ttB = [ttB; newttB];
c = [c newc];
time  = [time newtime];

%save('numsamples2.mat','x');
%save('numbonds2.mat','y') ;
%save('onrateconstants2.mat','z');

%close all;
%figure;
% curr2 = size(y,2);
% x =[1:curr2];
% subplot1 = subplot(2,2,1);
% plot(x,y(1,:));
% hold on;
% plot(x,y(2,:));
% plot(x,y(3,:));
% plot(x,y(4,:));
% plot(x,y(5,:));
% title('Times to Bond per sample');
% %legend1 = legend(subplot1,'show');
% %set(legend1,...
%     %'Position',[0.435366740294704 0.761163035904754 0.0947712408236689 0.137071647376658]);
% hold off;

%subplot2 = subplot(2,2,2);
for n = 1:4
    sumttB = 0;
    sumcount = 0;
    for i = 1:length(c)
        if ttB(i,n) < ttB(i,5)
            sumcount = sumcount + 1;
            sumttB = sumttB + ttB(i,n);
        end
    end
    kon(n) = 1/((sumttB/sumcount)*mean(c(:))); %angstroms^3/picosecond*atoms
    %hold on;
end
%plot([0:3], kon);
%title('On-rate constant vs. Initial # of bonds');
%hold off;

konMs = kon.*10^12*(6.02214*10^23)/(10^10)^3;
save('randwalk_output_test12.0.mat');


% x = [];
% for i = 1:11
%     c = 10*i-10;
%     [times,locs,types,rs,orients,vf1(i),tt,tD(i),eD(i),bP,nB,ttB] = RandWalk(18,18,5,2,c,10000,-1);
%     x(i) = i;
%     close all;
%     figure;
%     plot(vf1,tD);
%     hold on;
%     plot(vf1,eD);
% end

% x = [];
% 
% y = [];
% z = [];
% for n = 1:5
%     for i = 1:10
%         x(i) = i;
%         [times,locs,types,rs,orients,vf1,tt,tD,eD,bP,nB,ttB,c] = RandWalk(18,18,5,2,0,1000,n-1);
%         y(n,i) = ttB;
%         z(n,i) = c;
%     end
% end
% %close all;
% figure;
% subplot1 = subplot(2,2,1);
% plot(x,y(1,:));
% hold on;
% plot(x,y(2,:));
% plot(x,y(3,:));
% plot(x,y(4,:));
% plot(x,y(5,:));
% title('Times to Bond per sample');
% legend1 = legend(subplot1,'show');
% set(legend1,...
%     'Position',[0.435366740294704 0.761163035904754 0.0947712408236689 0.137071647376658]);
% hold off;
% 
% subplot2 = subplot(2,2,2);
% for n = 1:5
%     kon(n) = 1/(mean(y(n,:))*mean(z(n,:)));
%     hold on;
% end
% plot([0:4], kon);
% title('On-rate constant vs. Initial # of bonds');
% hold off;
% 
% subplot(2,2,3);
% hold on;
% bar([0,1,2,3,4],mean(y,2));
% errorbar([0,1,2,3,4],mean(y,2),std(y,0,2));
% title('Average Times to Bond by initial bonding');
% hold off;
% 
% subplot(2,2,4);
% plot([0,1,2,3,4],sum(z,2)/i);
% title('Average Number of Bonds by initial bonding');

% clear all;
% 
% x = [];
% for i = 1:5
%     x(i) = i;
%     [times,locs,types,rs,orients,vf1(i),tt,tD,eD,bP,nB,ttB] = RandWalk(18,18,5,2,50,2,-1);
% 
%     close all;
%     figure;
%     plot(x,vf1);
%     hold on;
% end