function [] = CPnl_net_figs()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for outputting Figures - Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

liqinj = 0;

if liqinj == 0
    load('IRF_CPnl_noLI.mat');
    load('varsCPnl_noLI.mat');
    LI = 'noliqinj';
elseif liqinj == 1
    load('IRF_CPnl_LI.mat');
    load('varsCPnl_LI.mat');
    LI = 'liqinj';
end

netstruct  = 'CPnl_';
fullsim    = 100;
numSC      = 4;

% Matrix preallocation

numvars   = size(var_list_,1); 
fields    = fieldnames(SC1);

finirfs_SC1  = zeros(numvars,fullsim);
finirfs_SC2  = zeros(numvars,fullsim);
finirfs_SC3  = zeros(numvars,fullsim);
finirfs_SC4  = zeros(numvars,fullsim);


if liqinj == 0   
    for i = 1:numvars
        finirfs_SC1(i,:)  = SC1.(fields{i})+SC1.(fields{i+numvars});
        finirfs_SC2(i,:)  = SC2.(fields{i})+SC2.(fields{i+numvars});
        finirfs_SC3(i,:)  = SC3.(fields{i})+SC3.(fields{i+numvars});
        finirfs_SC4(i,:)  = SC4.(fields{i})+SC4.(fields{i+numvars});
    end   
elseif liqinj == 1    
    for i = 1:numvars
        finirfs_SC1(i,:)  = SC1.(fields{i})+SC1.(fields{i+numvars})+SC1.(fields{i+(2*numvars)});
        finirfs_SC2(i,:)  = SC2.(fields{i})+SC2.(fields{i+numvars})+SC2.(fields{i+(2*numvars)});
        finirfs_SC3(i,:)  = SC3.(fields{i})+SC3.(fields{i+numvars})+SC3.(fields{i+(2*numvars)});
        finirfs_SC4(i,:)  = SC4.(fields{i})+SC4.(fields{i+numvars})+SC4.(fields{i+(2*numvars)});
    end
end

truncsim = 50;
tvec     = linspace(1,truncsim,truncsim);

finirf_mat(:,:,1) = finirfs_SC1(:,1:truncsim);
finirf_mat(:,:,2) = finirfs_SC2(:,1:truncsim);
finirf_mat(:,:,3) = finirfs_SC3(:,1:truncsim);
finirf_mat(:,:,4) = finirfs_SC4(:,1:truncsim);

trunc2sim = 10;   % Set if further truncating of the IRFs is needed (from visual inspection)

%----------------------------------------------------------------------------------------
% Figures for the paper
%----------------------------------------------------------------------------------------

fname_fig_S  = '/Users/nscholte/Desktop/Research/Ch.2 - MP Transmission/Presentations/Conference/Figures/CP'; % Slides: Directory for Figures;
fname_fig_P  = '/Users/nscholte/Desktop/Research/Ch.2 - MP Transmission/Drafts/v1/Figures/CP'; % Draft: Directory for Figures;


%% Interest rates

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank spreads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,2,1);
        h=plot(tvec(1:truncsim),finirf_mat(1,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(1,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(1,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(1,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}-r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','southoutside','Orientation','horizontal');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,2,2);
        h=plot(tvec(1:truncsim),finirf_mat(2,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(2,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(2,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(2,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AD}-r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [15 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 15 3]);
str1 = 'IBspreads_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));  

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prime lending rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(2,4,1);
         h=plot(tvec(1:truncsim),finirf_mat(3,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(3,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(3,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(3,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{l,A}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','southeast');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(4,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(4,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(4,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(4,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{l,B}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(5,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(5,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(5,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(5,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{l,C}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,4);
         h=plot(tvec(1:truncsim),finirf_mat(6,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(6,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(6,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(6,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{l,D}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Deposit rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(2,4,5);
         h=plot(tvec(1:truncsim),finirf_mat(7,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(7,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(7,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(7,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{d,A}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,6);
        h=plot(tvec(1:truncsim),finirf_mat(8,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(8,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(8,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(8,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{d,B}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,7);
        h=plot(tvec(1:truncsim),finirf_mat(9,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(9,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(9,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(9,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{d,C}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,8);
        h=plot(tvec(1:truncsim),finirf_mat(10,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(10,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(10,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(10,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{d,D}$'},'FontSize',8,'interpreter','latex');
        %legend({'$SC_{1}$','$SC_{2}$'},'FontSize',6,'Location','best','interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
        %legend('boxoff') 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 6]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 6]);
str1 = 'RErates_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));  

%% Interbank volumes

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank lending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,3,1);
        h=plot(tvec(1:truncsim),finirf_mat(11,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(11,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(11,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(11,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,BA}$'},'FontSize',8,'interpreter', 'latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,2);
        h=plot(tvec(1:truncsim),finirf_mat(12,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(12,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(12,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(12,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,AC}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,3);
        h=plot(tvec(1:truncsim),finirf_mat(13,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(13,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(13,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(13,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,AD}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
% subplot(2,3,4);
%         h=plot(tvec(1:truncsim),finirf_mat(14,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(14,1:truncsim,2)*100,'--b',...
%         tvec(1:truncsim),finirf_mat(14,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(14,1:truncsim,4)*100,'-.r');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,AB}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
% subplot(2,3,5);
%         h=plot(tvec(1:truncsim),finirf_mat(15,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(15,1:truncsim,2)*100,'--b',...
%         tvec(1:truncsim),finirf_mat(15,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(15,1:truncsim,4)*100,'-.r');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,CA}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
% subplot(2,3,6);
%         h=plot(tvec(1:truncsim),finirf_mat(16,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(16,1:truncsim,2)*100,'--b',...
%         tvec(1:truncsim),finirf_mat(16,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(16,1:truncsim,4)*100,'-.r');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,DA}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [15 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 15 3]);
str1 = 'IBvolumes_';
print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));    
  

%% Real Economy volumes
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firm credit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(2,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(17,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(17,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(17,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(17,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$X_{A}$'},'FontSize',8,'interpreter', 'latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(18,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(18,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(18,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(18,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$X_{B}$'},'FontSize',8,'interpreter', 'latex');        
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(19,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(19,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(19,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(19,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$X_{C}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(20,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(20,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(20,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(20,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$X_{D}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Household deposits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(2,4,5);
        h=plot(tvec(1:truncsim),finirf_mat(21,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(21,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(21,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(21,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$D_{A}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,6);
        h=plot(tvec(1:truncsim),finirf_mat(22,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(22,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(22,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(22,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$D_{B}$'},'FontSize',8,'interpreter', 'latex');        
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,7);
        h=plot(tvec(1:truncsim),finirf_mat(23,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(23,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(23,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(23,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$D_{C}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(2,4,8);
        h=plot(tvec(1:truncsim),finirf_mat(24,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(24,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(24,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(24,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$D_{D}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 6]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 6]);
str1 = 'REvolumes_';
print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));        
        
%% Defaults

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,3,1);
        h=plot(tvec(1:truncsim),finirf_mat(25,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(25,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(25,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(25,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,AB}$'},'FontSize',8,'interpreter', 'latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,2);
        h=plot(tvec(1:truncsim),finirf_mat(26,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(26,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(26,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(26,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,CA}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,3);
        h=plot(tvec(1:truncsim),finirf_mat(27,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(27,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(27,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(27,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,DA}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock');
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [15 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 15 3]);
str1 = 'IBdefaults_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

figure
subplot(1,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(28,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(28,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(28,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(28,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,A}$'},'FontSize',8,'interpreter','latex');
        legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(29,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(29,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(29,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(29,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,B}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(30,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(30,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(30,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(30,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,C}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(31,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(31,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(31,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(31,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,D}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'REdefaults_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));




%% Figures for appendix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Interbank rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure
subplot(1,3,1);
        h=plot(tvec(1:truncsim),finirf_mat(36,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(36,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(36,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(36,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','southeastoutside','Orientation','horizontal');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,2);
        h=plot(tvec(1:truncsim),finirf_mat(37,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(37,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(37,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(37,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,3,3);
        h=plot(tvec(1:truncsim),finirf_mat(38,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(38,1:truncsim,2)*100,'--b',...
        tvec(1:truncsim),finirf_mat(38,1:truncsim,3)*100,':.k',tvec(1:truncsim),finirf_mat(38,1:truncsim,4)*100,'-.r');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AD}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'IBrates_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Isolated banking shock %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure

for i = 1:3
subplot(3,3,(3*i)-2);
        h=plot(tvec(1:truncsim),finirf_mat(36,1:truncsim,i+1)*100-finirf_mat(36,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}','SC_{3}','SC_{4}'},'FontSize',6,'Location','southeastoutside','Orientation','horizontal');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(3,3,(3*i)-1);
        h=plot(tvec(1:truncsim),finirf_mat(37,1:truncsim,i+1)*100-finirf_mat(37,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(3,3,3*i);
        h=plot(tvec(1:truncsim),finirf_mat(38,1:truncsim,i+1)*100-finirf_mat(38,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AD}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 10]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 10]);
str1 = 'IBratesBSonly_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));

end