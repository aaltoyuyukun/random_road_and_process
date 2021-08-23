clear all;
clc;

w = 5;
max_h = 20;
amplitude = -max_h/2 + rand(w,1) * max_h;

T = 999;
P = (T+1)*10;
t = linspace(0,T,P); % 平均插值, linspace(0,999,10000) 999/10000 = 0.0999
t2 = t.';
len = numel(t);

ang_speed = rand(5,1);
phase = 2*pi*rand(5,1);

for j= 1:w
    for i= 1:len
        Wave_Load_Normal(j,i) = amplitude(j,:) * sin(ang_speed(j,:) * t2(i,:) + phase(j,:)); 
    end
end

Wave_Load_Random = sum(Wave_Load_Normal);

figure;
plot(t2,Wave_Load_Random);
xlabel("Time [s]");
ylabel("WaveHeight [m]");
title("Random Waves");
xlim([0 100]); % x axial from 0-100

wave_max = max(Wave_Load_Random);
wave_min = min(Wave_Load_Random);
gap = 1;
num_lev = (wave_max - wave_min) / gap;

for i = 1:num_lev
    layer(i) = wave_max - (i * gap);
    count(i) = sum(Wave_Load_Random > layer(i));
    cumula_prob(i) = count(i) / len;
    prob_vector(i) = wave_min + (i * gap);
end

num3 = numel(prob_vector) - 1;
k = zeros(1,num3);
for a = 1:num3
    k(a) = (cumula_prob(a+1) - cumula_prob(a))/(prob_vector(a+1)-prob_vector(a));
end
k(num3+1) = 0;

% scatter(prob_vector, cumula_prob);
% subplot(1,2,2);
% ployfit(prob_vector, cumula_prob);

% fit curve on calculated distribution points
x = 1:numel(k);
x2 = polyfit(prob_vector,k,9);
y2 = polyval(x2,prob_vector);
figure
scatter(prob_vector,k);
figure
plot(prob_vector,y2);

sam = T/len; %每个周期的时间
Fs = 1/sam;

[c,hist,edges,rmm,idx] = rainflow(Response(:,2),Fs);
T = array2table(c,'VariableNames',{'Count','Range','Mean','Start','End'});

figure;
histogram('BinEdges','edges','BinCounts',sum(hist,2));
title('Cycle counts as a function of stress');
xlabel('Stress Range');
ylabel('Cycle Counts');

figure;
rainflow(Response,Fs);

random_fft = fft(Wave_Load_Random);
% Compute the two-sided spectrum P2.
% Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length len.
p2 = abs(random_fft);
p1 = p2(1:len/2+1);
p1(2:end-1) = 2 * p1(2:end-1);

figure;
rainflow(Response,Fs);

% figure
% % % Define the frequency domain fre_dom_random
% fre_dom_random = Fs * (0:(len/2))/len;
% plot(fre_dom_random,p1);

% t_crit=1000;
% for i=1:w
%     for j=1:100
%         if j*2*pi()/ang_speed(i,:) < t_crit
%             num_waves_vector(j,:)=1;
%         else
%             num_waves_vector(j,:)=0;
%         end
%     end
%     number_fullwaves(i,:)=sum(num_waves_vector);
%     phase_crit(i,:)=pi()/2-((100-number_fullwaves(i,:)*2*pi()/ang_speed(i,:))*ang_speed(i,:));
%     for k=1:len
%         FreakLoadNormal(i,k)=amplitude(i,:)*sin(ang_speed(i,:)*t2(k,:)+phase_crit(i,:));
%         FreakLoadComponents(i,k)=amplitude(i,:)*sin(ang_speed(i,:)*t2(k,:)+phase_crit(i,:))+(i-1)*max_h;
%     end
%     freak_moment_value(i,:)=FreakLoadNormal(i,t_crit); 
%     if freak_moment_value(i,:)<0
%         for l=1:len
%             FreakLoadComponents(i,l)=(amplitude(i,:)*sin(ang_speed(i,:)*t2(l,:)+phase_crit(i,:)))*(-1)+(i-1)*max_h;
%             FreakLoadNormal(i,l)=amplitude(i,:)*sin(ang_speed(i,:)*t2(l,:)+phase_crit(i,:))*(-1);
%         end
%     end
% end
% 
