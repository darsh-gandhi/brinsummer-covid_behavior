% this code accompanies the SIR w/ Hospitalization + Death coupled group behavior
% model; the study focuses on the COVID outbreak across New York City from
% March 2nd 2020 (day_0) to Aug 27th (?) 2020 (day_179).

clear; close all;
colors = 1/255 * [56 182 255; 255 222 89; 0 191 99; 0 98 155; 189 151 0; 0 105 54; 255 49 49; 115 115 115];

%% data (post processing)

load('COVID_data_NYC_2020.mat', 'I_tot_NYC_daily_avg', 'D_tot_NYC_daily_avg', 'I_tot_NYC_daily', 'D_tot_NYC_daily', 'Dates');

hosp_cases = [1	2 7	2 14 8	8 18 37	60	78	106	152	177	210	338	385	494	631	738	734	773	1131 1301	1434	1535	1585	1394	1509	1858	1816	1768	1684	1712	1460	1406	1750	1590	1559	1353	1293	1048	951	1156	984	866	782	790	584	536	641	587	519	448	474	327	307	370	356	313	254	263	208	216	241	213	213	176	166	123	105	180	113	130	142	128	91	80	111	114	80	83	102	74	86	79	76	75	79	84	54	60	81	64	73	75	53	42	54	60	65	65	49	52	35	47	52	61	49	56	45	25	44	48	44	52	43	33	43	32	47	45	39	45	30	30	29	54	47	34	46	45	34	26	49	43	44	34	50	26	38	38	24	40	33	40	35	25	38	38	30	37	34	27	33	40	32	38	34	23	21	26	41	25	31	39	38	23	21	33	39	35	37	31	27	23	36	40	28	31	22	22	20	42	21	29	18	31	22	31	19	26	31	25	32	11	17	32	29	46	34	37	21	24	29	34	38	42	46	26	37	38	50	51	41	49	25	44	55	48	46	47	55	38	49	52	50	58	74	48	35	40	69	74	52	53	60	47	44	56	46	45	59	60	48	55	64	75	69	85	64	60	65	80	89	99	80	92	101	92	127	102	120	132	135	119	116	150	131	141	156	158	127	163	188	200	176	218	182	178	173	222	220	278	246	244	209	227	266	238	214	231	274	239	246	331	287	303	276	267	264	241	356	324	316	287];
timespan = 0:1:length(I_tot_NYC_daily_avg)-1; % days from 3/1 - 12/31/20; Mar 1 := Day 0, Dec 31 := Day 305

data_mat = [timespan; I_tot_NYC_daily_avg; hosp_cases; D_tot_NYC_daily_avg];
data_mat = data_mat(:,2:181);
data_mat(1,:) = data_mat(1,:)-1;
timespan = data_mat(1,:);

%% data visual

figure(1)
plot(data_mat(1,:), data_mat(2,:),'LineWidth',1.5,'Color',colors(1,:))
hold on
plot(data_mat(1,:),data_mat(3,:),'LineWidth',1.5,'Color',colors(2,:))
plot(data_mat(1,:),data_mat(4,:),'LineWidth',1.5,'Color',colors(3,:))
legend('Infections','Hospitalizations','Deaths')
xlim([0 180])

%% define params + states

prop_hosp_1 = sum(data_mat(3,:))/sum(data_mat(2,:)); % percent of infected people hospitalized

% parameter order:
% z := [(1) pi, (2) alpha, (3) r, (4) sigma, (5) p, (6) gamma, (7) delta,
% (8) xi, (9) psi_12^max, (10) psi_21^max, (11) a_1, (12) a_2, (13) beta]
z = [1/180, 0.05, prop_hosp_1, 0.15, 0.95, 0.1, 0.41, 0.05, 5, 5, 5e2, 5e2, 0.4];

% state order:
% y := [(1) S_1, (2) I_1, (3) R_1, (4) S_2, (5) I_2, (6) R_2, (7) H, (8) D];
y0 = [z(8)*8.8e6 z(8)*data_mat(2,1) 0 (1-z(8))*8.8e6 (1-z(8))*data_mat(2,1) 0 data_mat(3,1) data_mat(4,1)]; % assuming 5% of people comply before mandates introduced on march 12
N_0 = sum(y0);

% solve ode using MATLAB solver
[t,y] = ode15s(@(t,y) odesystem(t,y,z),timespan,y0);
% 
%% plot ode_sols

figure(2)
plot(t,y(:,1),'LineWidth',2.0,'Color',colors(1,:))
hold on
plot(t,y(:,2),'LineWidth',2.0,'Color',colors(2,:))
plot(t,y(:,3),'LineWidth',2.0,'Color',colors(3,:))
legend('S_1','I_1','R_1','Location','northeast')
xlabel('Time (days)')
ylabel('# of Individuals')
xlim([0 180])
hold off

figure(3)
plot(t,y(:,4),'LineWidth',2.0,'Color',colors(4,:))
hold on
plot(t,y(:,5),'LineWidth',2.0,'Color',colors(5,:))
plot(t,y(:,6),'LineWidth',2.0,'Color',colors(6,:))
legend('S_2','I_2','R_2','Location','east')
xlabel('Time (days)')
ylabel('# of Individuals')
xlim([0 180])
hold off

figure(4)
plot(t,y(:,7),'LineWidth',2.0,'Color',colors(7,:))
hold on
plot(t,y(:,8),'LineWidth',2.0,'Color',colors(8,:))
legend('H','D','Location','east')
xlabel('Time (days)')
ylabel('# of Individuals')
xlim([0 180])
hold off

%% data fitting

% parameter order:
% z := [(1) pi, (2) alpha, (3) r, (4) sigma, (5) p, (6) gamma, (7) delta,
% (8) xi, (9) psi_12^max, (10) psi_21^max, (11) a_1, (12) a_2, (13) beta]

lb = z;
ub = z;

lb(9:13) = [0,              0,         0,        0, 0];
ub(9:13) = [1e-2*N_0, 1e-2*N_0, 5e-1*N_0, 5e-1*N_0, 100];

%fitting only to infection case data:
data = data_mat(2,:);

omega = 1e2; nu=1;
[z(9), z(10), z(11), z(12), z(13)] = group_multistart(data,y0,timespan,lb,ub,omega,nu)
