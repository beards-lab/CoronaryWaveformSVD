%% OVERVIEW:

%Code and data associated with preprint:
%Shape Analysis of Coronary Flow Waveforms using Singular Value Decomposition
%Victoria E. Sturgess, Noah A. Schenk, Joshua W. Ziegele, Salman I. Essajee, 
% Johnathan D. Tune, Indika Rajapakse, C. Alberto Figueroa, Daniel A. Beard
%bioRxiv 2026.08.26.743980; doi: https://doi.org/10.64898/2026.08.26.743980


%% Load data
%Load data file
load('uploadData_30beats.mat')
Case="Baseline"; %always running with baseline case
%Pigs Included
IDs={'OSS1145', 'OSS1146','OSS1147','OSS1177','OSS1178','OSS1179','OSS1197','OSS1198','OSS1199','OSS1201','OSS1169','OSS1172','OSS1193','OSS1194','OSS1195','OSS1196','OSS1202','OSS1348','OSS1349','OSS1350', 'OSS1415', 'OSS1416','OSS1377','OSS1379','York1359', 'York1337', 'YorkD1', 'YorkD2', 'YorkD4', 'YorkD5', 'YorkD6', 'YorkD7'};
Types={'O','O','O','O','O','O','L','L','L','L','C','C','C','C','C','C','C','L','L','L','L','C','C','C','Y','Y','Y','Y','Y','Y','Y','Y'};
Types2=[3,3,3,3,3,3,2,2,2,2,1,1,1,1,1,1,1,2,2,2,2,1,1,1,5,4,4,4,4,4,4,4];

%Types:
% Obese Paced Ossabaw: O or 3
% Lean Paced Ossabaw: L or 2
% Lean Ossabaw (Control): C or 1
% Yorkshire (Control): Y 
    % Yorkshire Paced: 5
    % Yorkshire Non-paced: 4

%% Identify 30 beats, calculate CCV values
%%%%% Note-Section has been run already to create uploadedData %%%%%%

% CCV_value=zeros(1,length(IDs));
% CCV_per=zeros(1,length(IDs));
% for i=1:length(IDs)
%     PigNum=IDs{i};
%     Case_IDs=data.(PigNum).summary.Properties.RowNames;
%     for j=1:1 %only baseline for now
%         Case=convertCharsToStrings(Case_IDs(j));
%         Current=data.(PigNum).(Case);
% 
%         %find a location when dplv/dt is a maximum
%         dPLVdt=Current.PLV(2:end)-Current.PLV(1:end-1);
%         dPLVdt(isnan(dPLVdt)) = [];
%         dplv_dt_smooth = smoothdata(dPLVdt,"movmean",50);
%         mean_dPLVdt=mean(dplv_dt_smooth);
%         std_dPLVdt=std(dplv_dt_smooth);
%         peaks_dPLVdt= mean_dPLVdt+0.8*std_dPLVdt;
%         [pks,locs]=findpeaks(dplv_dt_smooth, 'MinPeakHeight',peaks_dPLVdt, 'MinPeakDistance',400);
% 
%         avg_peak_value=mean(pks); 
% 
%         locs2=locs;
%         if locs(1)<100
%             locs2(1)=[];
%             locs2=locs2(1:30);
%         end
% 
%         EDP_locs=zeros(length(locs2),1);
%         EDP_t=zeros(length(locs2),1);
%         EDP_values=zeros(length(locs2),1);
%         EDP_dPdt_values=zeros(length(locs2),1);
% 
%         for k=1:length(locs2)
%             [M,I]=min(abs(0.2*avg_peak_value-dplv_dt_smooth(locs2(k)-100:locs2(k))));
%             EDP_locs(k)=locs2(k)-101+I;
%             EDP_dPdt_values(k)=dplv_dt_smooth(locs2(k)-101+I);
%             EDP_t(k)=Current.t(locs2(k)-101+I)-Current.t(locs(1));
%             EDP_values(k)=Current.PLV(locs2(k)-101+I);
%         end
% 
%         t1=locs(1);
%         t2=locs(31)-1;
%         if mod(length(Current.t(t1:t2)),2)==1
%             t2=locs(31);
%         end
%         new_locs=locs(find(locs<=(t2+1)));
%         diff=new_locs(2:end)-new_locs(1:end-1);
%         HR=60/(0.001*mean(diff));
% 
%         CCV_value(i)=mean(abs(diff-mean(diff)));
%         CCV_per(i)=mean(abs(diff-mean(diff)))/mean(diff); 
% 
%         data_30beats.(PigNum).(Case).PLV=Current.PLV(t1:t2);
%         data_30beats.(PigNum).(Case).PAo=Current.PAo(t1:t2);
%         data_30beats.(PigNum).(Case).Q=Current.Q(t1:t2);
%         data_30beats.(PigNum).(Case).t=Current.t(t1:t2)-Current.t(t1);
%         data_30beats.(PigNum).(Case).HR=HR;
%         data_30beats.(PigNum).(Case).beats=30;
%         data_30beats.(PigNum).(Case).LVEDP=mean(EDP_values);
%     end
% end

%% Viewing Signal in the frequency domain

for i=1:length(IDs) %Looking at all pigs
% for i=24 %looking at individual pig
    PigNum=IDs{i};

    Current=data_30beats.(PigNum).(Case);
    t=Current.t;
    Q=Current.Q;    
    [~,F_Q,fshift]= calcFourierTransform(Q',t');

    F_Q2_m3=F_Q(fshift>=0);
    fshift2_m3=fshift(fshift>=0);
    mag3=abs(F_Q2_m3);
    ph3=angle(F_Q2_m3);

    data_30beats.(PigNum).(Case).raw_freq0=mag3(1);
    data_30beats.(PigNum).(Case).raw_mag=mag3;
    data_30beats.(PigNum).(Case).raw_ph=ph3;
    data_30beats.(PigNum).(Case).raw_real=real(F_Q2_m3);
    data_30beats.(PigNum).(Case).raw_im=imag(F_Q2_m3);
    data_30beats.(PigNum).(Case).raw_fshift=fshift2_m3;

    
    %Uncomment to create Figure 1 panel b
    % % Magnitude of complex number (real + imaginary component)
    % figure()
    % hold on
    % plot(fshift2_m3,mag3,'LineWidth',2,'color', 'k', 'LineStyle','-','DisplayName',"30 beats")
    % xlim([0 20])
    % title(PigNum)
    % xlabel('Frequency (Hz)');
    % ylabel('Magnitude');
    % box on
    % set(gca,'FontSize',18,'LineWidth',1.25)
    % % saveas(gcf,strcat("Results_Paper/FrequencyDomain_",PigNum,".png" ))
    % 
    % % Plotting only the real component
    % figure()
    % hold on
    % for i=1:15
    %     xline(fshift2_m3(i*30+1), '--')
    % end
    % plot(fshift2_m3,abs(real(F_Q2_m3)),'LineWidth',2,'color', 'r', 'LineStyle','-','DisplayName',"30 beats")
    % xlim([0 18])
    % title(PigNum)
    % xlabel('Frequency (Hz)');
    % ylabel('Real Component');
    % box on
    % set(gca,'FontSize',18,'LineWidth',1.25)
    % 
    % % Plotting only the imaginary component
    % figure()
    % hold on
    % for i=1:15
    %     xline(fshift2_m3(i*30+1), '--')
    % end
    % plot(fshift2_m3,abs(imag(F_Q2_m3)),'LineWidth',2,'color', 'r', 'LineStyle','-','DisplayName',"30 beats")
    % xlim([0 18])
    % title(PigNum)
    % xlabel('Frequency (Hz)');
    % ylabel('Imaginary Component');
    % box on
    % set(gca,'FontSize',18,'LineWidth',1.25)

end

%% Decide on number of harmonics
%%%%% Calculates SSE, MSE, and NMSE for the baseline signal %%%%%%
%%%%% for each pig using a waveform reconstructed with      %%%%%%
%%%%% 1,2,...,49,50 harmonics of the heart rate             %%%%%%

SSE_all=zeros(50,length(IDs));
MSE_all=zeros(50,length(IDs));
NMSE_all=zeros(50,length(IDs));
for i=1:length(IDs)
    % for i=1:1
    PigNum=IDs{i};
    Current=data_30beats.(PigNum).(Case);
    t=Current.t;
    Q=Current.Q;
    
    freq0=Current.raw_freq0;
    mag=Current.raw_mag;
    ph=Current.raw_ph;
    fshift=Current.raw_fshift;
    
    
    for N=1:50
        Q_recon=ones(1,length(t))*freq0;
        
        %loop for reconstructing waveform with N harmonics
        for s=1:N
            j=30*s+1;
            Q_recon=Q_recon+2*mag(j)*cos(2*pi*fshift(j)*t+ph(j));
        end
        SSE_all(N,i)=sum((Q-Q_recon').^2);
        MSE_all(N,i)=mean((Q-Q_recon').^2);
        NMSE_all(N,i)=mean((Q-Q_recon').^2)/std(Q);
        
    end
end

MSE_all_raw=MSE_all;
NMSE_all_raw=NMSE_all;


%Figure 2
figure('position', [100 100 900 350]);

subplot(1,2,1); hold on
mean_mse_FourierModes=mean(MSE_all_raw,2);
for i=1:length(IDs)
    plot(MSE_all_raw(:,i),'LineStyle','-')
end
plot(mean_mse_FourierModes,'k', 'LineWidth',4, 'LineStyle','-');

box on
xlim([0 50])
set(gca,'FontSize',15,'LineWidth',1.25)
ylabel("MSE")
xlabel("# of Fourier Modes")

subplot(1,2,2); hold on
mean_mse_FourierModes=mean(NMSE_all_raw,2);
for i=1:length(IDs)
    plot(NMSE_all_raw(:,i),'LineStyle','-')
end
plot(mean_mse_FourierModes,'k', 'LineWidth',4, 'LineStyle','-');

box on
xlim([0 50])
set(gca,'FontSize',15,'LineWidth',1.25)
ylabel("NMSE")
xlabel("# of Fourier Modes")

%% Calculate N=15 MSE
%%%%% Recalculate SSE, MSE, and NMSE using specifically N=15 %%%%%%

N=15;
SSE_15=zeros(1,length(IDs));
MSE_15=zeros(1,length(IDs));
NMSE_15=zeros(1,length(IDs));

for i=1:length(IDs)
% for i=2:2 %use to run for specific pigs if plotting figures
    PigNum=IDs{i};
    Current=data_30beats.(PigNum).(Case);
    t=Current.t;
    Q=Current.Q;

    freq0=Current.raw_freq0;
    mag=Current.raw_mag;
    ph=Current.raw_ph;
    fshift=Current.raw_fshift;

    Q_recon=ones(1,length(t))*freq0;

    %loop for reconstructing waveform with N harmonics
    for s=1:N
        j=30*s+1;
        Q_recon=Q_recon+2*mag(j)*cos(2*pi*fshift(j)*t+ph(j));
    end

    
    %Uncomment for panels A-C of Figure 3
    % figure(), hold on
    % plot(t,Q,'k', 'LineWidth',3, 'DisplayName',"Raw Signal")
    % plot(t,Q_recon,'r','LineWidth',3,'LineStyle',':', 'DisplayName',"Reconstructed Signal")
    % % title(strcat(int2str(N)," Harmonics"))
    % title(PigNum)
    % xlabel('Time (seconds)');
    % ylabel('Flow (ml/min)');
    % legend('Location','southeast')
    % xlim([0 10*60/Current.HR])
    % set(gca,'FontSize',18,'LineWidth',1.25)
    % box on
    
    % MSE_cur=mean(abs(Q-Q_recon'))
    SSE_15(i)=sum((Q-Q_recon').^2);
    MSE_15(i)=mean((Q-Q_recon').^2);
    NMSE_15(i)=mean((Q-Q_recon').^2)/std(Q);

end

% Figure 3 pannels d and e
figure('Position',[200,50,350,700]);
subplot(2,1,1);hold on
scatter(CCV_value(Types2==1), NMSE_15(Types2==1),100,'red','filled','MarkerEdgeColor', 'k');
scatter(CCV_value(Types2==2), NMSE_15(Types2==2),100,'green','filled','MarkerEdgeColor', 'k');
scatter(CCV_value(Types2==3), NMSE_15(Types2==3),100,'blue','filled','MarkerEdgeColor', 'k');
scatter(CCV_value(Types2==4), NMSE_15(Types2==4),100,'red',"^",'filled','MarkerEdgeColor', 'k');
scatter(CCV_value(Types2==5), NMSE_15(Types2==5),100,'green',"^",'filled','MarkerEdgeColor', 'k');
ylabel("NMSE")
xlabel("CCDV (ms)")
set(gca,'FontSize',15,'LineWidth',1.25)
box on
xlim([0,25])
subplot(2,1,2);hold on
scatter(100*CCV_per(Types2==1), NMSE_15(Types2==1),100,'red','filled','MarkerEdgeColor', 'k');
scatter(100*CCV_per(Types2==2), NMSE_15(Types2==2),100,'green','filled','MarkerEdgeColor', 'k');
scatter(100*CCV_per(Types2==3), NMSE_15(Types2==3),100,'blue','filled','MarkerEdgeColor', 'k');
scatter(100*CCV_per(Types2==4), NMSE_15(Types2==4),100,'red',"^",'filled','MarkerEdgeColor', 'k');
scatter(100*CCV_per(Types2==5), NMSE_15(Types2==5),100,'green',"^",'filled','MarkerEdgeColor', 'k');
ylabel("NMSE")
xlabel("CCDV (% of cardiac cycle)")
set(gca,'FontSize',15,'LineWidth',1.25)
box on

% Linear regression fits vor NMSE vs CCV value and percent
disp("Linear regression fit for CCV_value vs NMSE_15")
fitlm(CCV_value,NMSE_15, 'poly1')

disp("Linear regression fit for CCV_per vs NMSE_15")
fitlm(CCV_per,NMSE_15, 'poly1')

%% Create Matrix A 
%%%%% size 31x Npigs                                   %%%%%%
%%%%% contains average flow (zero frequency component) %%%%%%
%%%%% contains average flow (zero frequency component) %%%%%%
%%%%% and magnitude and phase from first 15 harmonics  %%%%%%

N_harmonics=15;
matrixA_raw=zeros(N_harmonics*2,length(IDs));

for i=1:length(IDs)
    PigNum=IDs{i};
    Current=data_30beats.(PigNum).(Case);
    
    matrixA_raw(1,i)=Current.raw_freq0;
    
    idx=1;
    for s=1:N_harmonics
        j=30*s+1;
        matrixA_raw(1+idx,i)=Current.raw_real(j);
        matrixA_raw(1+idx+N_harmonics,i)=Current.raw_im(j);

        idx=idx+1;
    end

    
end

%% Create average signal
% computed assuming HR of 60 beats per minute 
% or a cardiac cycle frequency of 1 hz
avg_feat=mean(matrixA_raw,2);
feat=complex(avg_feat(1:15),avg_feat(16:30));
mag=abs(feat);
ph=angle(feat);
fshift=1; %1 hz
t=0:0.01:5;
Q_recon=zeros(1,length(t));
for j=1:15
    Q_recon=Q_recon+2*mag(j)*cos(2*pi*j*fshift*t+ph(j));
end

% figure()
% plot(t,Q_recon,"Color",'k','Linewidth',1)
% set(gca,'FontSize',18,'LineWidth',1.25)


%% SVD and scree plot 
%%%%% performs svd calculation and creates scree plot and %%%%%
%%%%% a percent variance explained plot %%%%%%
[leftSV_raw,Sigma_raw,rightSV_raw]=svd(matrixA_raw);
singularvalues=zeros(1,length(IDs));
for i=1:(min([1+2*N_harmonics,length(IDs)]))
    singularvalues(i)=Sigma_raw(i,i);
end
singularvalues_raw=singularvalues;



figure()
plot(1:length(singularvalues),singularvalues,"-o","Color",'k','Linewidth',2)
set(gca,'FontSize',18,'LineWidth',1.25)
xlabel("SVD mode")
ylabel("Singular Value")

var_ex=zeros(1,length(IDs));
total_var=sum(singularvalues.^2);
for i=1:length(IDs)
    var_ex(i)=sum(singularvalues(1:i).^2)/total_var;
end

% Figure 4 left panel
figure(), hold on
plot(1:length(var_ex),100*var_ex,"-o","Color",'k','Linewidth',2)
scatter(5,100*var_ex(5),100,'r','filled')
set(gca,'FontSize',18,'LineWidth',1.25)
box on
xlabel("Number of SVD modes")
ylabel("Percent Variance Explained")

%% Coefficients and modal contributions
%%%%% plots for the mean and average contirbutions of each mode %%%%%%
%%%%% modal contributions are calculated as    %%%%%%
%%%%% singular values * right singular vectors %%%%%%

Coeff=Sigma_raw*transpose(rightSV_raw);

Coeff_mean=mean(Coeff');
Coeff_stdev=std(Coeff');

% Figure 4 right panel
figure()
errorbar(Coeff_mean,Coeff_stdev,'k','LineWidth',1.25)
set(gca,'FontSize',18,'LineWidth',1.25)
box on
xlabel("SVD Mode")
ylabel("Modal Contributions")
ylim([-70,10])


%% Plotting effects of each SVD mode of shape of the coronary waveform
%%%%% compares how +/- 1 stdev of each SVD mode impacts the %%%%%%
%%%%% coronary flow waveform compared to the average coronary waveform %%%%%%

Coeff=Sigma_raw*transpose(rightSV_raw);

Coeff_stdev=std(Coeff');

for sv=1:5
    avg_feat=mean(matrixA_raw,2);
    feat_plus_sv=avg_feat+1*Coeff_stdev(sv)*leftSV_raw(:,sv);
    feat_minus_sv=avg_feat-1*Coeff_stdev(sv)*leftSV_raw(:,sv);

    C1=complex(avg_feat(2:16),avg_feat(17:31));
    C2=complex(feat_plus_sv(2:16),feat_plus_sv(17:31));
    C3=complex(feat_minus_sv(2:16),feat_minus_sv(17:31));


    mag1=abs(C1);
    ph1=angle(C1);
    mag2=abs(C2);
    ph2=angle(C2);
    mag3=abs(C3);
    ph3=angle(C3);
    
    fshift=1; %1 hz
    t=0:0.01:3;
    Q_recon1=ones(1,length(t))*avg_feat(1);
    Q_recon2=ones(1,length(t))*feat_plus_sv(1);
    Q_recon3=ones(1,length(t))*feat_minus_sv(1);
    
    
    for j=1:N_harmonics
        Q_recon1=Q_recon1+2*mag1(j)*cos(2*pi*j*fshift*t+ph1(j));
        Q_recon2=Q_recon2+2*mag2(j)*cos(2*pi*j*fshift*t+ph2(j));
        Q_recon3=Q_recon3+2*mag3(j)*cos(2*pi*j*fshift*t+ph3(j));
    end

    %Figure 5 left panels
    figure(), hold on
    plot(t,Q_recon1,"Color",'k','Linewidth',2.5)
    plot(t,Q_recon2,"Color",'r','Linewidth',2.5, 'LineStyle',':')
    plot(t,Q_recon3,"Color",'b','Linewidth',2.5, 'LineStyle',':')
    set(gca,'FontSize',18,'LineWidth',1.25)
    box on
    xlabel("Cardiac Cycles")
    ylabel("Flow (ml/min)")
    ylim([-10 110])
    % ylim([-1.5 1.5])
    xlim([0 2])
    title(strcat("Mode ", num2str(sv)))
end

%% Load Table with Experimental Data
BaselinePig=readtable('uploadTable.csv',"ReadRowNames",true);

% uncomment to run if any changes to the selected 30 beats
% for i=1:length(IDs)
%     PigNum=IDs{i};
%     PLV=data_30beats.(PigNum).(Case).PLV;
%     dPLVdt=(PLV(2:end)-PLV(1:end-1))/0.001;
% 
%     BaselinePig{PigNum, "HR"}=data_30beats.(PigNum).(Case).HR;
%     BaselinePig{PigNum, "MBP"}=mean(data_30beats.(PigNum).(Case).PAo);
%     BaselinePig{PigNum, "PLVmax"}=max(data_30beats.(PigNum).(Case).PLV);
%     BaselinePig{PigNum, "PLVmin"}=min(data_30beats.(PigNum).(Case).PLV);
%     BaselinePig{PigNum, "LVEDP"}=data_30beats.(PigNum).(Case).LVEDP;
%     BaselinePig{PigNum, "AvgFlow"}=mean(data_30beats.(PigNum).(Case).Q);
%     BaselinePig{PigNum, "AvgFlow_g"}=BaselinePig{PigNum, "AvgFlow"}/BaselinePig{PigNum, "Heart_Weight"};
%     BaselinePig{PigNum, "FlowPulse"}=max(data_30beats.(PigNum).(Case).Q)-min(data_30beats.(PigNum).(Case).Q);
%     BaselinePig{PigNum, "dPLVdtmax"}=max(dPLVdt);
%     BaselinePig{PigNum, "dPLVdtmin"}=min(dPLVdt);
% 
%     BaselinePig{PigNum, "CCVval"}=CCV_value(i);
%     BaselinePig{PigNum, "CCVper"}=CCV_per(i);
%     BaselinePig{PigNum, "MSE"}=MSE_15(i);
% end

%% Scatterplots of SVD and Experimental Data
%%%%% used for visual comparison between pig-specific component in a  %%%%%%
%%%%% right singular vector and experimental data %%%%%%
%%%%% alter modes and col_names to plot other comparisons  %%%%%%

numPigs=length(IDs);
F=zeros(1,sum(numPigs));
RSV_zvalue_raw=(rightSV_raw-mean(rightSV_raw))./std(rightSV_raw);

%Use for figure 5 right side
modes=[1,2,3,4,5];
col_names=["AvgFlow","HeartRate","MBP","LVEDP","AvgFlow"];

%Use for figure 7
% modes=[1,1,3,3];
% col_names=["AvgFlow","MaxHyperemicFlow","AvgFlow","MaxHyperemicFlow",];

%Use for figure S1
% modes=[4];
% col_names=["Endo_EpiRatio"];

for i=1:length(modes)
    for j=1:length(IDs)
        F(j)=BaselinePig.(col_names(i))(j);
    end
    figure();hold on
    scatter(F(Types2==1), RSV_zvalue_raw(Types2==1,modes(i)),100,'red','filled','MarkerEdgeColor', 'k');
    scatter(F(Types2==2), RSV_zvalue_raw(Types2==2,modes(i)),100,'green','filled','MarkerEdgeColor', 'k');
    scatter(F(Types2==3), RSV_zvalue_raw(Types2==3,modes(i)),100,'blue','filled','MarkerEdgeColor', 'k');
    scatter(F(Types2==4), RSV_zvalue_raw(Types2==4,modes(i)),100,'red',"^",'filled','MarkerEdgeColor', 'k');
    scatter(F(Types2==5), RSV_zvalue_raw(Types2==5,modes(i)),100,'green',"^",'filled','MarkerEdgeColor', 'k');
    % title(strcat("Singular Vector ",num2str(i)))
    ylabel(strcat("Mode ",num2str(modes(i)), ": z-value"))
    xlabel(col_names(i))
    set(gca,'FontSize',18,'LineWidth',1.25)
    box on
end

%% Correlations between singular vectors and experimental data
%%%%% runs a linear regression correlation analysis between  %%%%%%
%%%%% each pair of experimental data and right singular vectors (modes 1:5) %%%%%%

sv = 5; %Look at top 5 SVD modes
correlationAnalysis = struct();

paramList = fieldnames(BaselinePig);
paramListSkips = {'Properties', 'Row', 'Variables'};
validParamList = ~ismember(paramList,paramListSkips);
paramList = paramList(validParamList);

% Loop through each parameter
for s = 1:length(paramList)

    %loop through each singular value
    numSingularValues = min(sv, length(Sigma_raw));
    results = zeros(numSingularValues, 5);
    for i = 1:numSingularValues
        % extracts right singualr value
        y  = RSV_zvalue_raw(:,i);
        %extracts param
        x = BaselinePig.(paramList{s})(:);

        %remove NaN
        validIdx = ~(isnan(x) | isnan(y));
        x = x(validIdx);
        y = y(validIdx);
        
        %normalize data
        xmin = min(x);
        xmax = max(x);
        x = (x - xmin) ./ (xmax - xmin);

        %skip if too little points 
        if numel(x) < 6
            results(i,:) = [NaN, NaN, NaN, NaN, NaN];
            continue;
        end

        %Linear fit
        p = polyfit(x, y, 1);
        y_fit = polyval(p, x);
        [r, pValue] = corr(x, y, 'Type', 'Pearson');
        
        % R² calculation
        SS_res = sum((y - y_fit).^2);
        SS_tot = sum((y - mean(y)).^2);
        R2 = 1 - (SS_res / SS_tot);

        % Store results in a row
        results(i, :) = [p(1), p(2), R2, r, pValue];
    end
    %create table for the results 
    rowNames = arrayfun(@(i) sprintf('SingularValue%d', i), 1:numSingularValues, 'UniformOutput', false);
    tbl = array2table(results, 'VariableNames', {'Slope', 'Intercept', 'R2', 'PearsonCorrelation', 'PValue'}, ...
        'RowNames', rowNames);
    %store results 
    correlationAnalysis.(paramList{s}) = tbl;
end

%% Create bubble plot (figure 6)
%%%%% visual representation of p-value and PearsonCorrelation %%%%%%
%%%%% from the analysis in the previous section %%%%%%

%changing names to fit matlab formatting
paramList2=paramList;
paramList2{9,1} = 'HW-normAvgFlow';
paramList2{14,1} = 'Endo-EpiRatio';
paramList2{17,1} = '\DeltaHyperemicFlow';

%Store p-value and PearsonCorrelation for numSingularValues SVD modes and
%each experimental parameter
numSummaryParams = 20;
corrMatrix_p = nan(numSingularValues, numSummaryParams);
corrMatrix_r = nan(numSingularValues, numSummaryParams);

matvalue_p = 'PValue';
matvalue_r = 'PearsonCorrelation';

for s = 1:numSummaryParams
    
    % Retrieve the table of regression results for the current summary parameter
    tbl = correlationAnalysis.(paramList{s});

    % For each singular value (up to maximumNumRSV), store the Pearson correlation
    for i = 1:numSingularValues
        corrMatrix_p(i, s) = tbl.(matvalue_p)(i);
        corrMatrix_r(i, s) = tbl.(matvalue_r)(i);
    end
end

%Create figure 8
figure('Position', [100, 100, 1000 500]);
mainAx = axes('Position', [0.1 0.45 0.8 0.5]);  % Leave space for legend on the right

[x, y] = meshgrid(1:numSummaryParams,1:numSingularValues);

dotSize= -log10(corrMatrix_p(:));
dotSize(dotSize>20)=20;
dotColor = corrMatrix_r(:);

%Create colormap
n = 256; % number of colors in colormap
half = floor(n/2);

blue2white = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
white2red = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
customCMap = [blue2white; white2red];

scatter(x(:),y(:),dotSize*30,dotColor,'filled')

colormap(customCMap); 
caxis([-1 1]);        
colorbar;
% grid on;

set(gca,'XTick', 1:numSummaryParams, 'XTickLabel', paramList2);
set(gca,'YDir','reverse', 'YTick', 1:numSingularValues, 'YTickLabel', [1,2,3,4,5],'FontSize',15);
ylabel("SVD Mode",'FontSize',18)
xlim([0.5,numSummaryParams+0.5])
ylim([0.5,5.5])
box on


% Create invisible axes for dot size legend
dotAx = axes('Position', [0.9 0.5 0.1 0.5]);  % Adjust position as needed
axis(dotAx, 'off');  % Hide axis lines and ticks
hold(dotAx, 'on');

dotSizes = [20, 10, 5, 2, 1];
yLoc = linspace(0.8, 0.2, numel(dotSizes));
xLoc = 0.5;  % Center of the legend axis


for i = 1:numel(dotSizes)
    scatter(dotAx, xLoc, yLoc(i), dotSizes(i)*30, [0.5 0.5 0.5], ...
            'filled', 'MarkerEdgeColor', 'k');
end


%% Construct representative waveforms for the three groups of ossabaw
%Figure 8- 
%%%% Creates a waveform representing each experimental ossabaw group using %%%%%
%%%% the average waveform and group-average contributions of specific SVDmodes %%%%%
%%%% Alter include_modes to change which modal contributions are %%%%%
%%%% represented by the specific-group average                   %%%%%

Coeff_OssControl=Coeff(:,Types2==1);
Coeff_OssLeanPaced=Coeff(:,Types2==2);
Coeff_OssObesePaced=Coeff(:,Types2==3);

%each group starts off being assigned modal coefficients for average pig
Coeff_mean=mean(Coeff,2);
avgCoeff_OssControl=mean(Coeff_OssControl,2);
avgCoeff_OssLeanPaced=mean(Coeff_OssLeanPaced,2);
avgCoeff_OssObesePaced=mean(Coeff_OssObesePaced,2);

%these modes are updated to represent the average over a specific group
%(lean control, lean paced, obese paced)

include_modes=[1,2] %for select modes
% include_modes=1:15 %for all modes

curCoeff_OssControl=Coeff_mean;
curCoeff_OssLeanPaced=Coeff_mean;
curCoeff_OssObesePaced=Coeff_mean;

% update contributions of include_modes to represent averages from lean
% control, lean paced, and obese paced
for idx=1:length(include_modes)
    sv_i=include_modes(idx);
    curCoeff_OssControl(sv_i)=avgCoeff_OssControl(sv_i);
    curCoeff_OssLeanPaced(sv_i)=avgCoeff_OssLeanPaced(sv_i);
    curCoeff_OssObesePaced(sv_i)=curCoeff_OssObesePaced(sv_i);
end

%multiple by the left singular values to get a vector for each pig which
%can be used for waveform reconstruction
feat_OssControl=leftSV_raw*curCoeff_OssControl;
feat_OssLeanPaced=leftSV_raw*curCoeff_OssLeanPaced;
feat_OssObesePaced=leftSV_raw*curCoeff_OssObesePaced;

%waveform reconstruction
C1=complex(feat_OssControl(2:16),feat_OssControl(17:31));
C2=complex(feat_OssLeanPaced(2:16),feat_OssLeanPaced(17:31));
C3=complex(feat_OssObesePaced(2:16),feat_OssObesePaced(17:31));

mag1=abs(C1);
ph1=angle(C1);
mag2=abs(C2);
ph2=angle(C2);
mag3=abs(C3);
ph3=angle(C3);

fshift=1; %1 hz
t=0:0.01:3;

Q_recon1=ones(1,length(t))*feat_OssControl(1);
Q_recon2=ones(1,length(t))*feat_OssLeanPaced(1);
Q_recon3=ones(1,length(t))*feat_OssObesePaced(1);

for j=1:N_harmonics
    Q_recon1=Q_recon1+2*mag1(j)*cos(2*pi*j*fshift*t+ph1(j));
    Q_recon2=Q_recon2+2*mag2(j)*cos(2*pi*j*fshift*t+ph2(j));
    Q_recon3=Q_recon3+2*mag3(j)*cos(2*pi*j*fshift*t+ph3(j));
end

figure(), hold on
plot(t,Q_recon1,"Color",'r','Linewidth',1.5)
plot(t,Q_recon2,"Color",'g','Linewidth',1.5)
plot(t,Q_recon3,"Color",'b','Linewidth',1.5)
set(gca,'FontSize',18,'LineWidth',1.25)
box on
xlabel("Cardiac Cycles")
ylabel("Flow (ml/min)")
ylim([0 90])
% ylim([-1.5 1.5])
xlim([0 2])

%% Stats test on RSV z-values between Ossabaw groups
%%%% Performs a t-test to compare contributions from a right -singularvector %%%%%
%%%% between different Ossawbaw groups %%%%%
%%%% Performed for SVD modes 1-5 %%%%%
%%%% Part of table 1 in manuscript %%%%%


RSV_zvalue_raw=(rightSV_raw-mean(rightSV_raw))./std(rightSV_raw);

RSV_zvalue_OssControl=RSV_zvalue_raw(Types2==1,:);
RSV_zvalue_OssLeanPaced=RSV_zvalue_raw(Types2==2,:);
RSV_zvalue_OssObesePaced=RSV_zvalue_raw(Types2==3,:);

avgRSV_zvalue_OssControl=mean(RSV_zvalue_OssControl,1);
avgRSV_zvalue_OssLeanPaced=mean(RSV_zvalue_OssLeanPaced,1);
avgRSV_zvalue_OssObesePaced=mean(RSV_zvalue_OssObesePaced,1);

stdevRSVz_OssControl=std(RSV_zvalue_OssControl,1);
stdevRSVz_OssLeanPaced=std(RSV_zvalue_OssLeanPaced,1);
stdevRSVz_OssObesePaced=std(RSV_zvalue_OssObesePaced,1);

for SVD_mode_i=1:5
    fprintf("For SVD mode %i \n",SVD_mode_i)

    fprintf("The group mean and stdev for Ossabaw Lean Control is %f and %f \n", mean(RSV_zvalue_OssControl(:,SVD_mode_i)), std(RSV_zvalue_OssControl(:,SVD_mode_i)))
    fprintf("The group mean and stdev for Ossabaw Lean Paced is %f and %f \n", mean(RSV_zvalue_OssLeanPaced(:,SVD_mode_i)), std(RSV_zvalue_OssLeanPaced(:,SVD_mode_i)))
    fprintf("The group mean and stdev for Ossabaw Obese Paced is %f and %f \n", mean(RSV_zvalue_OssObesePaced(:,SVD_mode_i)), std(RSV_zvalue_OssObesePaced(:,SVD_mode_i)))
    
    
    [h,p_CL,ci,stats]=ttest2(RSV_zvalue_OssControl(:,SVD_mode_i),RSV_zvalue_OssLeanPaced(:,SVD_mode_i));
    fprintf("t-test between Ossabaw Lean Control and Ossabaw Lean Paced gives a p-value of %f \n",p_CL)
    
    [h,p_CO,ci,stats]=ttest2(RSV_zvalue_OssControl(:,SVD_mode_i),RSV_zvalue_OssObesePaced(:,SVD_mode_i));
    fprintf("t-test between Ossabaw Lean Control and Ossabaw Obese Paced gives a p-value of %f \n",p_CO)
    
    [h,p_LO,ci,stats]=ttest2(RSV_zvalue_OssLeanPaced(:,SVD_mode_i),RSV_zvalue_OssObesePaced(:,SVD_mode_i));
    fprintf("t-test between Ossabaw Lean Paced and Ossabaw Obese Paced gives a p-value of %f \n",p_LO)

end
%% Stats test on average experimental data between Ossabaw groups
%%%% Performs a t-test to compare select experiment data %%%%%
%%%% between different Ossawbaw groups %%%%%
%%%% Part of table 1 in manuscript %%%%%

col_names=["AvgFlow","HeartRate","MBP","LVEDP","AvgFlow"];
numPigs=length(IDs);
F=zeros(1,sum(numPigs));

for i=1:length(col_names)
    for j=1:length(IDs)
        F(j)=BaselinePig.(col_names(i))(j);
    end
    fprintf("For experimental column %s \n",col_names(i))

    fprintf("The group mean and stdev for Ossabaw Lean Control is %f and %f \n", mean(F(Types2==1)), std(F(Types2==1)))
    fprintf("The group mean and stdev for Ossabaw Lean Paced is %f and %f \n", mean(F(Types2==2)), std(F(Types2==2)))
    fprintf("The group mean and stdev for Ossabaw Obese Paced is %f and %f \n", mean(F(Types2==3)), std(F(Types2==3)))
    
    [h,p_CL,ci,stats]=ttest2(F(Types2==1),F(Types2==2));
    fprintf("\nt-test between Ossabaw Lean Control and Ossabaw Lean Paced gives a p-value of %f \n",p_CL)
    
    [h,p_CO,ci,stats]=ttest2(F(Types2==1),F(Types2==3));
    fprintf("t-test between Ossabaw Lean Control and Ossabaw Obese Paced gives a p-value of %f \n",p_CO)
    
    [h,p_LO,ci,stats]=ttest2(F(Types2==2),F(Types2==3));
    fprintf("t-test between Ossabaw Lean Paced and Ossabaw Obese Paced gives a p-value of %f \n",p_LO)
end


%% Predict hyperemic flow
%%%%% Comparisons between true hyperemic flow measurements and   %%%%%%
%%%%% hyperemic flow measurements predicted by linear regression %%%%%%
%%%%% to average baseline flow, SVD mode 1, SVD mode 3           %%%%%%
%%%%% and a multiple linear regression to SVD modes 1 and 3      %%%%%%

 CFRs_true=zeros(1,sum(numPigs));
 Hyperemic_true=zeros(1,sum(numPigs));
 BaselineFlow_true=zeros(1,sum(numPigs));
 MBP_true=zeros(1,sum(numPigs));

 for i=1:length(IDs)
    CFRs_true(i)=BaselinePig.("CFR")(i);
    Hyperemic_true(i)=BaselinePig.("MaxHyperemicFlow")(i);
    BaselineFlow_true(i)=BaselinePig.("AvgFlow")(i);
    MBP_true(i)=BaselinePig.("MBP")(i);
end

 CFRs_true2=CFRs_true(~isnan(Hyperemic_true));
 Hyperemic_true2=Hyperemic_true(~isnan(Hyperemic_true));
 BaselineFlow_true2=BaselineFlow_true(~isnan(Hyperemic_true));
 MBP_true2=MBP_true(~isnan(Hyperemic_true));

 raw_rsv1=rightSV_raw(~isnan(Hyperemic_true),1);
 raw_rsv3=rightSV_raw(~isnan(Hyperemic_true),3);

 p1 = polyfit(BaselineFlow_true2, Hyperemic_true2, 1);
 Hyperemic_from_baseline = polyval(p1, BaselineFlow_true2);
 [r1, pValue1] = corr(BaselineFlow_true2, Hyperemic_true2, 'Type', 'Pearson');
 SS_res1 = sum((Hyperemic_true2 - Hyperemic_from_baseline).^2);
 SS_tot = sum((Hyperemic_true2 - mean(Hyperemic_true2)).^2);
 R2_1 = 1 - (SS_res1 / SS_tot);
 MSE1=SS_res1/length(Hyperemic_true2);
 fprintf("Linear regression between Baseline Flow and Hyperemic Flow has a R^2 value of %f and MSE of %f \n",R2_1,MSE1)

 p2 = polyfit(raw_rsv1, Hyperemic_true2, 1);
 Hyperemic_from_rsv1 = polyval(p2, raw_rsv1);
 [r2, pValue2] = corr(raw_rsv1, Hyperemic_true2', 'Type', 'Pearson');
 SS_res2 = sum((Hyperemic_true2 - Hyperemic_from_rsv1').^2);
 SS_tot = sum((Hyperemic_true2 - mean(Hyperemic_true2)).^2);
 R2_2 = 1 - (SS_res2 / SS_tot);
 MSE2=SS_res2/length(Hyperemic_true2);
 fprintf("Linear regression between RSV1 and Hyperemic Flow has a R^2 value of %f and MSE of %f \n",R2_2,MSE2)

 p3 = polyfit(raw_rsv3, Hyperemic_true2, 1);
 Hyperemic_from_rsv3 = polyval(p3, raw_rsv3);
 [r3, pValue3] = corr(raw_rsv3, Hyperemic_true2', 'Type', 'Pearson');
 SS_res3 = sum((Hyperemic_true2 - Hyperemic_from_rsv3').^2);
 SS_tot = sum((Hyperemic_true2 - mean(Hyperemic_true2)).^2);
 R2_3 = 1 - (SS_res3 / SS_tot);
 MSE3=SS_res3/length(Hyperemic_true2);
 fprintf("Linear regression between RSV3 and Hyperemic Flow has a R^2 value of %f and MSE of %f \n",R2_3,MSE3)

 %fit with curve fitting app in matlab
 Hyperemic_from_rsv1rsv3 = 48.14 + -312.6*raw_rsv1 + 141.4*raw_rsv3;
 SS_res4 = sum((Hyperemic_true2 - Hyperemic_from_rsv1rsv3').^2);
 SS_tot = sum((Hyperemic_true2 - mean(Hyperemic_true2)).^2);
 R2_4 = 1 - (SS_res4 / SS_tot);
 MSE4=SS_res4/length(Hyperemic_true2);
 fprintf("Linear regression between RSV1&3 and Hyperemic Flow has a R^2 value of %f and MSE of %f \n",R2_4,MSE4)

 figure()
 subplot(2,2,1), hold on
 scatter(Hyperemic_from_baseline,Hyperemic_true2,40,'k', 'filled')
 plot(0:1:200,0:1:200, 'k', 'LineStyle','--')
 xlabel("Predicted Hyperemic Flow (mL/min)")
 ylabel("True Hyperemic Flow (mL/min)")
 set(gca,"FontSize",12)

 subplot(2,2,2), hold on
 scatter(Hyperemic_from_rsv1,Hyperemic_true2,40,'k', 'filled')
 plot(0:1:200,0:1:200, 'k', 'LineStyle','--')
 xlabel("Predicted Hyperemic Flow (mL/min)")
 ylabel("True Hyperemic Flow (mL/min)")
 set(gca,"FontSize",12)

 subplot(2,2,3), hold on
 scatter(Hyperemic_from_rsv3,Hyperemic_true2,40,'k', 'filled')
 plot(0:1:200,0:1:200, 'k', 'LineStyle','--')
 xlabel("Predicted Hyperemic Flow (mL/min)")
 ylabel("True Hyperemic Flow (mL/min)")
 set(gca,"FontSize",12)

 subplot(2,2,4), hold on
 scatter(Hyperemic_from_rsv1rsv3,Hyperemic_true2,40,'k', 'filled')
 plot(0:1:200,0:1:200, 'k', 'LineStyle','--')
 xlabel("Predicted Hyperemic Flow (mL/min)")
 ylabel("True Hyperemic Flow (mL/min)")
 set(gca,"FontSize",12)


%% Functions
% Calculate Fourier Modes
function [Fn,Fn2,fshift]= calcFourierTransform(f,t)
    a=size(t);
    N=a(1);
    if N==1
        disp("Size or orientation of t may be wrong, try transpose")
    end
    T=t(N);
    w=2*pi/T;
    fs=N/T;

    Fn=fft(f(1:a(1)))/N; %consider the last point (not identical to the first
    Fn2 = fftshift(Fn);
    fshift = (-N/2:N/2-1)*(fs/N);
    f_space=fshift(2)-fshift(1);
end

