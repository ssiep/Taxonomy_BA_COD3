%
%

clear all
clear classes
close all
%
% figure
% hold on

%%----------Initial conditions for equilibrium continuation---------------%
%-------------------------------------------------------------------------%

% Find analytical solution for equilibria----------------------------------
% (they depend only on two cartesian parameters: mu1, mu2)

syms vx vmu1 vmu2 vrr vph vth

%Zeroes_of_X_cartesian=solve(-vx^3+vmu2*vx+vmu1,vx)
vmu1=-vrr*sin(vth)*sin(vph);    % go to spherical coordinates
vmu2= vrr*sin(vth)*cos(vph);

% find solutions for x (y will always be equal to 0 at equilibria)
Zeroes_of_X=solve(-vx^3+vmu2*vx+vmu1,vx) %the first two sol. are the foci, the third is the saddle

% As initial values for the spherical parameters we will use th=3, ph=0,
% while rr will be allowed to vary to build each spherical shell
th_in=3; ph_in=0;
%rr_in=0.5:0.1:1;
rr_in=0.5;
rr_n=length(rr_in);

Zeroes_of_X_rr=subs(Zeroes_of_X,[vth, vph],[3,0]);

%create vector of initial conditions for varying values of rr only
InitCond_for_EquilCont=vpa(subs(Zeroes_of_X_rr,vrr,rr_in)); %remember to take only real part of existent solution (tiny imaginary part due to numerical computation)




%%-------Start Matcont Command line---------------------------------------%
init

bb_in=1; % this parameter of the model will not be changed in the script

for rr_index=1:rr_n
    
    %%-------Equilibria Continuation--------------------------------------%
    % follow the first focus-----------------------------------------------
    
    p=[rr_in(rr_index); th_in; ph_in; bb_in];
    ap=2;
    x0_in=eval(real(InitCond_for_EquilCont(1,rr_index)));
    
    [x0, v0]=init_EP_EP(@Cod3FocusTimeFor_SphericalParam,[x0_in; 0],p,ap);
    opt=contset;
    opt=contset(opt,'Singularities',1);
    [x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
    
    %%-------HOPF curve---------------------------------------------------%
    % follow the first H point found in EP_EP------------------------------
    
    ap=[2 3]; % activate parameters th and ph
    
    % find coordinates and parameters of first H point and initialize H_H
    line_H=find(strcmp({s.label},'H '));
    index_firstH=s(line_H(1)).index;
    
    p=[rr_in(rr_index); x(3,index_firstH); ph_in; bb_in];
    [x0,v0] = init_H_H(@Cod3FocusTimeFor_SphericalParam, [x(1,index_firstH); x(2,index_firstH)], p, ap)
    
    % continue the curve
    opt=contset;
    opt=contset(opt,'Singularities',1);
    opt=contset(opt,'MaxNUmPoints',500);
    opt=contset(opt,'TestTolerance',1e-5 );
    opt=contset(opt,'VarTolerance',1e-7 );
    [x_H,v_H,s_H,h_H,f_H]=cont(@hopf,x0,[],opt);
    
    
    line_BT=find(strcmp({s_H.label},'BT'));
    index_secBT=s_H(line_BT(2)).index;
    index_thiBT=s_H(line_BT(3)).index;
    hopf_curve{rr_index,1} = x_H(:,index_secBT:index_thiBT);
    s_hopf_curve{rr_index,1}= s_H;
    
    % transform to cartesian coordinates and save - removing the part
    % inside the fold betwen the two TB which is wrong
    
    hopf_curve{rr_index,2}(1:2,:) = hopf_curve{rr_index,1}(1:2,:); % x and y coordinates
    hopf_curve{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), hopf_curve{rr_index,1}(3,:), hopf_curve{rr_index,1}(4,:)); % nu,mu1,mu2
    
    % plot
    %         plot3(hopf_curve{rr_index,2}(3,:), hopf_curve{rr_index,2}(4,:), hopf_curve{rr_index,2}(5,:),'g')
    %
    %%------FOLD curve----------------------------------------------------%
    % follow second BT point found in H_H----------------------------------
    
    % find coordinates and parameters of second BT point and initialize
    % BT_LP
    line_BT=find(strcmp({s_H.label},'BT'));
    index_secBT=s_H(line_BT(2)).index;
    
    p=[rr_in(rr_index); x_H(3,index_secBT); x_H(4,index_secBT); bb_in];
    [x0,v0] = init_BT_LP(@Cod3FocusTimeFor_SphericalParam, [x_H(1,index_secBT); x_H(2,index_secBT)], p, ap)
    
    % continue the curve
    opt=contset;
    opt=contset(opt,'Singularities',1);
    [x_LP,v_LP,s_LP,h_LP,f_LP]=cont(@limitpoint,x0,v0,opt);
    
    % save
    fold{rr_index,1} = x_LP;
    s_fold{rr_index,1}= s_LP;
    
    % transform to cartesian coordinates and save
    fold{rr_index,2}(1:2,:) = fold{rr_index,1}(1:2,:); % x and y coordinates
    fold{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), fold{rr_index,1}(3,:), fold{rr_index,1}(4,:)); % nu,mu1,mu2
    
    %         % plot
    %         plot3(fold{rr_index,2}(3,:), fold{rr_index,2}(4,:), fold{rr_index,2}(5,:),'Color', [1 128/255 0])
    %
    %
    %%------FOLD of CYCLES curves-----------------------------------------%
    % follow all the GH points in H_H--------------------------------------
    
    % find rows of GH points
    line_GH=find(strcmp({s_H.label},'GH'));
    GH_n=size(line_GH,2);
    GH_n_min=1;
    GH_n_max=3;
    
    for i=1
        i
        %initialize GH_LPC
        p=[rr_in(rr_index); x_H(3,s_H(line_GH(i)).index); x_H(4,s_H(line_GH(i)).index); bb_in];
        [x0,v0] = init_GH_LPC(@Cod3FocusTimeFor_SphericalParam, x_H, p,  s_H(line_GH(i)),ap, 20,4,0.01)
        
        % continue the curve
        opt=contset;
        opt=contset(opt,'Singularities',1);
        opt=contset(opt,'MaxStepsize',1.0);
        %opt=contset(opt,'MaxNumPoints',800);
        [x_LPCf,v_LPCf,s_LPCf,h_LPCf,f_LPCf]=cont(@limitpointcycle,x0,v0,opt);
        
        % save
        fold_cycle{rr_index,1,i} = x_LPCf;
        s_fold_cycle{rr_index,1,i}= s_LPCf;
        
        index_LPC_th = size(x_LPCf,1)-1;
        index_LPC_ph = size(x_LPCf,1);
        
        % transform to cartesian coordinates and save
        fold_cycle{rr_index,2,i}(1:2,:) = fold_cycle{rr_index,1,i}(1:2,:); % x and y coordinates
        fold_cycle{rr_index,2,i}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), fold_cycle{rr_index,1,i}(index_LPC_th,:), fold_cycle{rr_index,1,i}(index_LPC_ph,:)); % nu,mu1,mu2
        
        %         % backward
        %         opt=contset;
        %         opt=contset(opt,'Singularities',1);
        %         opt=contset(opt,'Backward',1);
        %         opt=contset(opt,'MaxStepsize',2.0);
        %         opt=contset(opt,'MaxNumPoints',300);
        %         [x_LPCb,v_LPCb,s_LPCb,h_LPCb,f_LPCb]=cont(@limitpointcycle,x0,[],opt);
        %
        %         % save
        %         fold_cycle_b{rr_index,1,i} = x_LPCb;
        %         s_fold_cycle_b{rr_index,1,i}= s_LPCb;
        %
        %         index_LPC_th = size(x_LPCb,1)-1;
        %         index_LPC_ph = size(x_LPCb,1);
        %
        %         % transform to cartesian coordinates and save
        %         fold_cycle_b{rr_index,2,i}(1:2,:) = fold_cycle_b{rr_index,1,i}(1:2,:); % x and y coordinates
        %         fold_cycle_b{rr_index,2,i}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), fold_cycle_b{rr_index,1,i}(index_LPC_th,:), fold_cycle_b{rr_index,1,i}(index_LPC_ph,:)); % nu,mu1,mu2
        
        %plot
        %                 hold on
        %                 plot3(fold_cycle{rr_index,2,i}(3,:), fold_cycle{rr_index,2,i}(4,:), fold_cycle{rr_index,2,i}(5,:),'b')
        %
        %             if(i == 4)
        %                save('test')
        %end
        %end
    end
    %
    %     clear all
    %     load('test')
    %
    
    %
    %
    %
    
    %----------HOMOCLINIC from BT1----------------------------------------%
    
    index_firBT=s_H(line_BT(1)).index;
    
    p=[rr_in(rr_index); x_H(3,index_firBT); x_H(4,index_firBT); bb_in];
    %   [x,v] = init_BT_Hom(odefile, x, s, p, ap, ntst, ncol,TTolerance ,amplitude, extravec)
    [x0,v0] = init_BT_Hom(@Cod3FocusTimeFor_SphericalParam, [x_H(1,index_firBT); x_H(2,index_firBT)], s_H(line_BT(1)), p, ap, 80,4, 1e-5, 0.01, [0 1 1])
    
    
    % continue the curve
    opt=contset;
    %  opt=contset(opt,'Backward',1);
    %  opt=contset(opt,'Singularities',1);
    [x_HomTB1,v_HomTB1,s_HomTB1,h_HomTB1,f_HomTB1]=cont(@homoclinic,x0,[],opt);
    
    % save
    HomTB1{rr_index,1} = x_HomTB1;
    %s_HomTB1{rr_index,1}= s_HomTB1;
    
    % transform to cartesian coordinates and save
    HomTB1{rr_index,2}(1:2,:) = HomTB1{rr_index,1}(1:2,:); % x and y coordinates
    HomTB1{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), HomTB1{rr_index,1}(645,:), HomTB1{rr_index,1}(646,:)); % nu,mu1,mu2
    
    %         % plot
    %         hold on
    %         plot3(HomTB1{rr_index,2}(3,:), HomTB1{rr_index,2}(4,:), HomTB1{rr_index,2}(5,:),'c')
    %
    %----------HOMOCLINIC from BT2----------------------------------------%
    
    index_secBT=s_H(line_BT(2)).index;
    
    p=[rr_in(rr_index); x_H(3,index_secBT); x_H(4,index_secBT); bb_in];
    %   [x,v] = init_BT_Hom(odefile, x, s, p, ap, ntst, ncol,TTolerance ,amplitude, extravec)
    [x0,v0] = init_BT_Hom(@Cod3FocusTimeFor_SphericalParam, [x_H(1,index_secBT); x_H(2,index_secBT)], s_H(line_BT(1)), p, ap, 80,4, 1e-5, 0.02, [0 1 1]);
    
    
    % continue the curve
    opt=contset;
    %  opt=contset(opt,'Backward',1);
    %  opt=contset(opt,'Singularities',1);
    [x_HomTB2,v_HomTB2,s_HomTB2,h_HomTB2,f_HomTB2]=cont(@homoclinic,x0,[],opt);
    
    % save
    HomTB2{rr_index,1} = x_HomTB2;
    %s_HomTB2{rr_index,1}= s_HomTB2;
    
    % transform to cartesian coordinates and save
    HomTB2{rr_index,2}(1:2,:) = HomTB2{rr_index,1}(1:2,:); % x and y coordinates
    HomTB2{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), HomTB2{rr_index,1}(645,:), HomTB2{rr_index,1}(646,:)); % nu,mu1,mu2
    
    %         % plot
    %         hold on
    %         plot3(HomTB2{rr_index,2}(3,:), HomTB2{rr_index,2}(4,:), HomTB2{rr_index,2}(5,:),'c')
    %
    
    %----------HOMOCLINIC MIDDLE curve------------------------------------%
    
    index_GH=1;
    index_lastLPC_from4thGH=s_fold_cycle{rr_index,1,index_GH}(end).index;
    
    p=[rr_in(rr_index); fold_cycle{rr_index,1,index_GH}(164,index_lastLPC_from4thGH); fold_cycle{rr_index,1,index_GH}(165,index_lastLPC_from4thGH); bb_in];
    % [x,v] = init_LC_Hom(odefile, x, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)
    [x0,v0] = init_LC_Hom(@Cod3FocusTimeFor_SphericalParam, fold_cycle{rr_index,1,index_GH}, s_fold_cycle{rr_index,1,index_GH}(end), p, ap, 80, 4,[0 1 1],fold_cycle{rr_index,1,index_GH}(163,index_lastLPC_from4thGH)/2,0.01,0.01);
    
    % continue the curve
    opt=contset;
    %  opt=contset(opt,'Backward',1);
    %  opt=contset(opt,'Singularities',1);
    opt=contset(opt,'Adapt',1);
    opt=contset(opt,'MaxNumPoints',200);
    
    [x_HomMiddle,v_HomMiddle,s_HomMiddle,h_HomMiddle,f_HomMiddle]=cont(@homoclinic,x0,[],opt);
    
    % save
    HomMiddle{rr_index,1} = x_HomMiddle;
    %s_HomMiddle{rr_index,1}= s_HomMiddle;
    
    % transform to cartesian coordinates and save
    HomMiddle{rr_index,2}(1:2,:) = HomMiddle{rr_index,1}(1:2,:); % x and y coordinates
    HomMiddle{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), HomMiddle{rr_index,1}(645,:), HomMiddle{rr_index,1}(646,:)); % nu,mu1,mu2
    
    %
    %        % plot
    %        hold on
    %         plot3(HomMiddle{rr_index,2}(3,:), HomMiddle{rr_index,2}(4,:), HomMiddle{rr_index,2}(5,:),'c')
    %
    
    % backward
    
    p=[rr_in(rr_index); fold_cycle{rr_index,1,index_GH}(164,index_lastLPC_from4thGH); fold_cycle{rr_index,1,index_GH}(165,index_lastLPC_from4thGH); bb_in];
    % [x,v] = init_LC_Hom(odefile, x, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)
    [x0,v0] = init_LC_Hom(@Cod3FocusTimeFor_SphericalParam, fold_cycle{rr_index,1,index_GH}, s_fold_cycle{rr_index,1,index_GH}(end), p, ap, 80, 4,[0 1 1],fold_cycle{rr_index,1,index_GH}(163,index_lastLPC_from4thGH)/2,0.01,0.01);
    
    % continue the curve
    opt=contset;
    opt=contset(opt,'Backward',1);
    opt=contset(opt,'Adapt',1);
    opt=contset(opt,'MaxNumPoints',200);
    %  opt=contset(opt,'Singularities',1);
    [x_HomMiddleb,v_HomMiddleb,s_HomMiddleb,h_HomMiddleb,f_HomMiddleb]=cont(@homoclinic,x0,[],opt);
    
    % save
    HomMiddleb{rr_index,1} = x_HomMiddleb;
    %s_HomMiddleb{rr_index,1}= s_HomMiddleb;
    
    % transform to cartesian coordinates and save
    HomMiddleb{rr_index,2}(1:2,:) = HomMiddleb{rr_index,1}(1:2,:); % x and y coordinates
    HomMiddleb{rr_index,2}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), HomMiddleb{rr_index,1}(645,:), HomMiddleb{rr_index,1}(646,:)); % nu,mu1,mu2
    
    %         % plot
    %         hold on
    %         plot3(HomMiddleb{rr_index,2}(3,:), HomMiddleb{rr_index,2}(4,:), HomMiddleb{rr_index,2}(5,:),'c')
    %
    
    for i=2:3
        i
        %initialize GH_LPC
        p=[rr_in(rr_index); x_H(3,s_H(line_GH(i)).index); x_H(4,s_H(line_GH(i)).index); bb_in];
        [x0,v0] = init_GH_LPC(@Cod3FocusTimeFor_SphericalParam, x_H, p,  s_H(line_GH(i)),ap, 20,4,0.01)
        
        % continue the curve
        opt=contset;
        opt=contset(opt,'Singularities',1);
        opt=contset(opt,'MaxStepsize',1.0);
        %opt=contset(opt,'MaxNumPoints',800);
        [x_LPCf,v_LPCf,s_LPCf,h_LPCf,f_LPCf]=cont(@limitpointcycle,x0,v0,opt);
        
        % save
        fold_cycle{rr_index,1,i} = x_LPCf;
        s_fold_cycle{rr_index,1,i}= s_LPCf;
        
        index_LPC_th = size(x_LPCf,1)-1;
        index_LPC_ph = size(x_LPCf,1);
        
        % transform to cartesian coordinates and save
        fold_cycle{rr_index,2,i}(1:2,:) = fold_cycle{rr_index,1,i}(1:2,:); % x and y coordinates
        fold_cycle{rr_index,2,i}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), fold_cycle{rr_index,1,i}(index_LPC_th,:), fold_cycle{rr_index,1,i}(index_LPC_ph,:)); % nu,mu1,mu2
        
        %         % backward
        %         opt=contset;
        %         opt=contset(opt,'Singularities',1);
        %         opt=contset(opt,'Backward',1);
        %         opt=contset(opt,'MaxStepsize',2.0);
        %         opt=contset(opt,'MaxNumPoints',300);
        %         [x_LPCb,v_LPCb,s_LPCb,h_LPCb,f_LPCb]=cont(@limitpointcycle,x0,[],opt);
        %
        %         % save
        %         fold_cycle_b{rr_index,1,i} = x_LPCb;
        %         s_fold_cycle_b{rr_index,1,i}= s_LPCb;
        %
        %         index_LPC_th = size(x_LPCb,1)-1;
        %         index_LPC_ph = size(x_LPCb,1);
        %
        %         % transform to cartesian coordinates and save
        %         fold_cycle_b{rr_index,2,i}(1:2,:) = fold_cycle_b{rr_index,1,i}(1:2,:); % x and y coordinates
        %         fold_cycle_b{rr_index,2,i}(3:5,:) = focusSpericalToCartesian(rr_in(rr_index), fold_cycle_b{rr_index,1,i}(index_LPC_th,:), fold_cycle_b{rr_index,1,i}(index_LPC_ph,:)); % nu,mu1,mu2
        
        %plot
        %                 hold on
        %                 plot3(fold_cycle{rr_index,2,i}(3,:), fold_cycle{rr_index,2,i}(4,:), fold_cycle{rr_index,2,i}(5,:),'b')
        %
        %             if(i == 4)
        %                save('test')
        %end
        %end
    end
    
    %%------Close rr loop-----------------------------------------------------%
end

%%------ PLOT CURVES -----------------------------------------------------%

figure
% try
for index_plot = 1:rr_n
    % plot bif sets in parameter space
    hold on
    plot3(fold{index_plot,2}(3,:), fold{index_plot,2}(4,:), fold{index_plot,2}(5,:),'Color', [1 128/255 0])
    plot3(hopf_curve{index_plot,2}(3,:), hopf_curve{index_plot,2}(4,:), hopf_curve{index_plot,2}(5,:),'g')
    plot3(hopf_curve{index_plot,2}(3,:), hopf_curve{index_plot,2}(4,:), hopf_curve{index_plot,2}(5,:),'g')
    
    plot3(HomTB1{index_plot,2}(3,:), HomTB1{index_plot,2}(4,:), HomTB1{index_plot,2}(5,:),'c')
    for i=1:3
        plot3(fold_cycle{index_plot,2,i}(3,:), fold_cycle{index_plot,2,i}(4,:), fold_cycle{index_plot,2,i}(5,:),'b')
        % plot3(fold_cycle_b{index_plot,2,i}(3,:), fold_cycle_b{index_plot,2,i}(4,:), fold_cycle_b{index_plot,2,i}(5,:),'b')
    end
    plot3(HomTB2{index_plot,2}(3,:), HomTB2{index_plot,2}(4,:), HomTB2{index_plot,2}(5,:),'c')
    nmax=170;
    plot3(HomMiddle{index_plot,2}(3,:), HomMiddle{index_plot,2}(4,:), HomMiddle{index_plot,2}(5,:),'c')
    
    plot3(HomMiddleb{index_plot,2}(3,:), HomMiddleb{index_plot,2}(4,:), HomMiddleb{index_plot,2}(5,:),'c')
    
    title('Cod 3 Focus type')
    xlabel('mu1')
    ylabel('mu2')
    zlabel('nu')
    legend('Fold','Hopf','Homoclinic to saddle','Fold of Cycles')
end
% end


filename = datestr(now);
save = (filename);



%hold off