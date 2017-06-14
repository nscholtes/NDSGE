function [] = cycl_net_figs()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for outputting Figures - Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

liqinj = 0;

if liqinj == 0
    load('IRF_cyclnet_noLI.mat');
    load('varscyclnet_noLI.mat');
    LI = 'noliqinj';
elseif liqinj == 1
    load('IRF_cyclnet_LI.mat');
    load('varscyclnet_LI.mat');
    LI = 'liqinj';
end

fullsim = 100;
netstruct  = 'Cycl_';
numSC      = 2;

% Matrix preallocation

numvars   = size(var_list_,1); 
fields    = fieldnames(SC1);

finirfs_SC1  = zeros(numvars,fullsim);
finirfs_SC2  = zeros(numvars,fullsim);

% Use of certainty equivalence principle to obtain total IRF (= sum of individual IRFs)

if liqinj == 0
    
    for i = 1:numvars
        finirfs_SC1(i,:)  = SC1.(fields{i})+SC1.(fields{i+numvars});
        finirfs_SC2(i,:)  = SC2.(fields{i})+SC2.(fields{i+numvars});
    end
    
elseif liqinj == 1   
    for i = 1:numvars
        finirfs_SC1(i,:)  = SC1.(fields{i})+SC1.(fields{i+numvars})+SC1.(fields{i+(2*numvars)});
        finirfs_SC2(i,:)  = SC2.(fields{i})+SC2.(fields{i+numvars})+SC2.(fields{i+(2*numvars)});
    end
end

truncsim = 50;
tvec     = linspace(1,truncsim,truncsim);

finirf_mat(:,:,1) = finirfs_SC1(:,1:truncsim);
finirf_mat(:,:,2) = finirfs_SC2(:,1:truncsim);

trunc2sim = 10;   % Set if further truncating of the IRFs is needed (from visual inspection)

%----------------------------------------------------------------------------------------
% Figures for the paper
%----------------------------------------------------------------------------------------

fname_fig_S  = '/Users/nscholte/Desktop/Research/Ch.2 - MP Transmission/Presentations/Conference/Figures/Cycl'; % Slides: Directory for Figures;
fname_fig_P  = '/Users/nscholte/Desktop/Research/Ch.2 - MP Transmission/Drafts/v1/Figures/Cycl'; % Draft: Directory for Figures;

%% Interest rates

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank rate spreads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(1,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(1,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}-r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %l=legend({'Agg. productivity shock','+banking shock (A)'},'FontSize',6,'Location',...
            %'southoutside','interpreter','latex','Orientation','horizontal');
        %set(l,'units','pixels');
        %set(l,'position',[5 5 300 15])
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(2,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(2,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,BA}-r^{b,DB}$'},'FontSize',6,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(3,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(3,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,CD}-r^{b,AC}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(4,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(4,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,DB}-r^{b,CD}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','southoutside','Orientation','horizontal');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'IBspreads_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI))); 

%% Interbank volumes

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank lending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(13,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(13,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,AC}$'},'FontSize',8,'interpreter', 'latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(14,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(14,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,BA}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(15,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(15,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,CD}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(16,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(16,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$L^{b,DB}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
% subplot(2,4,5);
%         h=plot(tvec(1:truncsim),finirf_mat(17,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(17,1:truncsim,2)*100,'--b');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,AB}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
% subplot(2,4,6);
%         h=plot(tvec(1:truncsim),finirf_mat(18,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(18,1:truncsim,2)*100,'--b');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,BD}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
% subplot(2,4,7);
%         h=plot(tvec(1:truncsim),finirf_mat(19,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(19,1:truncsim,2)*100,'--b');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,CA}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
% subplot(2,4,8);
%         h=plot(tvec(1:truncsim),finirf_mat(20,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(20,1:truncsim,2)*100,'--b');
%         %grid on
%         set(gca,'FontSize',6,'box','off');
%         set(h,'linewidth',1.5);
%         title({'$B^{b,DC}$'},'FontSize',8,'interpreter', 'latex');
%         %xlabel('quarters after shock')
%         xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'IBvolumes_';
print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));    
  
        
%% Defaults

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
subplot(1,4,1);
        h=plot(tvec(1:truncsim),(finirf_mat(29,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(29,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,AB}$'},'FontSize',8,'interpreter', 'latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','northeast');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),(finirf_mat(30,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(30,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,BD}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),(finirf_mat(31,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(31,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,CA}$'},'FontSize',8,'interpreter', 'latex');
        %xlabel('quarters after shock');
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),(finirf_mat(32,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(32,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\delta^{b,DC}$'},'FontSize',6,'interpreter', 'latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'IBrepay_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
figure
subplot(1,4,1);
        h=plot(tvec(1:truncsim),(finirf_mat(33,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(33,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,A}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','best');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),(finirf_mat(34,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(34,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,B}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),(finirf_mat(35,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(35,1:truncsim,2))*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$\alpha^{f,C}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),(finirf_mat(36,1:truncsim,1))*100,'-k',tvec(1:truncsim),(finirf_mat(36,1:truncsim,2))*100,'--b');
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
str1 = 'RErepay_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI)));
        
%% Figures for appendix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Interbank rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure
subplot(1,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(41,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(41,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','southeast');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(42,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(42,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,BA}$'},'FontSize',6,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(43,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(43,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,CD}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(44,1:truncsim,1)*100,'-k',tvec(1:truncsim),finirf_mat(44,1:truncsim,2)*100,'--b');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,DB}$'},'FontSize',8,'interpreter','latex');
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
subplot(1,4,1);
        h=plot(tvec(1:truncsim),finirf_mat(41,1:truncsim,2)*100-finirf_mat(41,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,AC}$'},'FontSize',8,'interpreter','latex');
        %legend({'SC_{1}','SC_{2}'},'FontSize',6,'Location','northeast');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,2);
        h=plot(tvec(1:truncsim),finirf_mat(42,1:truncsim,2)*100-finirf_mat(42,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,BA}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,3);
        h=plot(tvec(1:truncsim),finirf_mat(43,1:truncsim,2)*100-finirf_mat(43,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,CD}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
subplot(1,4,4);
        h=plot(tvec(1:truncsim),finirf_mat(44,1:truncsim,2)*100-finirf_mat(44,1:truncsim,1)*100,'-k');
        %grid on
        set(gca,'FontSize',6,'box','off');
        set(h,'linewidth',1.5);
        title({'$r^{b,DB}$'},'FontSize',8,'interpreter','latex');
        %xlabel('quarters after shock')
        xlim([0,truncsim]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 3]);
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 3]);
str1 = 'IBratesBSonly_';

print(gcf,'-depsc2',fullfile(fname_fig_S,strcat(str1,netstruct,LI)));
print(gcf,'-depsc2',fullfile(fname_fig_P,strcat(str1,netstruct,LI))); 

close all;

end