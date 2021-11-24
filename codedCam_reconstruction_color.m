%% Coded Illumination for improved lensless imaging
% This is a demo code for reproducing some of the main results in
% "Coded illumination for improved lensless imaging"

clear; 
% close all;

% ---- SCENE OPTION ----

% the id for testing scene:
% 1 for K card, 
% 2 for Cookie box
result_id = 1;

% ----------------

% illum_type = 'random'; % we have 3 illumination types to test: random, uniform, shifting_dots
% numIllum = 49; % number of illumination patterns applied


fig1 = figure(1); clf;
set(fig1,'Renderer', 'painters', 'Position', [600 300 800 500]);
fig2 = figure(2); clf;
set(fig2,'Renderer', 'painters', 'Position', [100 300 500 400]);

cnt = 1;
for illum_type = ["uniform","random","shifting dots"]
    switch illum_type
        case 'uniform'
            illum_ind = 3; numIllum = 1;
        case 'random'
            illum_ind = 1; numIllum = 49;
        case 'shifting dots'
            illum_ind = 4; numIllum = 49;
    end

    result_folder = fullfile('data',sprintf('exp%03d',result_id));
    result_path = fullfile(result_folder,sprintf('data_coded_color%03d_nillum%03d_object',illum_ind,numIllum));
    
    % read data of the captured measurements with illumination patterns and the patterns themselves
    load(result_path); 
    
    % read data of calibrated system matrix
    load(fullfile(result_folder ,'system_matrix.mat')); 

    data_coded = im2double(data_coded);
    measBg = im2double(measBg);

    % want to view captured images? use the command below
    % implay(data_coded/max(data_coded(:)))

    %% Reconstruction with illumination patterns
    % variables names:
    % phi_left, phi_right: left and right system matrices
    % PL, PR: separable illumination patterns in matrix
    % data_coded: captured lensless measurements under multiple different illuminations

    % data_coded = im2double(w);
    [M,N] = size(phi_left);
    xhat_color = zeros(N,N,3);
    if illum_ind == 3
        lambda_reg = 5e-2;
    else
        lambda_reg = 1e-3;
    end
    for cc = 1:3
        xhat_color(:,:,cc) = mat2gray((recon_op_illum(phi_left, phi_right, PL,PR, permute(squeeze(data_coded(:,:,cc,:)) - measBg(:,:,cc),[3,1,2]) ,lambda_reg)));
    end

    % display reconstruction results for each pattern type
    figure(1);
    subplot(1,3,cnt);
    imagesc(xhat_color);
    colormap gray;
    axis image off;
    title(sprintf('%02d %s pattern(s)', numIllum, illum_type),'FontSize',14);
    % imwrite(xhat_color, sprintf('img%02d_pattern%02d_num%02d.png',  result_id, illum_ind,numIllum));
    cnt = cnt + 1;

    % Compare singular values of the left system matrices
    AL = (phi_left'*phi_left).*(PL*PL');
    SL = svd(AL); 
    SL = SL/SL(1);
    figure(2); hold on
    plot(SL,'LineWidth',2)
    title('Singular values of the (left) system matrix','FontSize',16)
end
legend('1 uniform','49 random', '49 shifting dots','FontSize',16)
%% functions for reconstruction
function x = recon_op_illum(phi_L,phi_R,PL,PR,y,lambda)
[M,N] = size(phi_L);
ncolor = size(y,4);

AL = (phi_L'*phi_L).*(PL*PL');
AR = (phi_R'*phi_R).*(PR*PR');

[UL,SL,VL] = svd(AL);
[UR,SR,VR] = svd(AR);
sL = diag(SL);
sR = diag(SR);

nillumL = size(PL,2);
nillumR = size(PR,2);

Yh = zeros(N,N,ncolor);
for ii = 1:nillumL
    for jj = 1:nillumR
        for cc = 1:ncolor
            Yh(:,:,cc) = Yh(:,:,cc) + diag(PL(:,ii))*phi_L'*squeeze(y(sub2ind([nillumL,nillumR],ii,jj),:,:,cc))*phi_R*diag(PR(:,jj));
        end
    end
end

x = zeros(N,N,ncolor);
for cc = 1:ncolor
    x(:,:,cc) = VL*((VL'*Yh(:,:,cc)*VR) ./ (sL*sR' + lambda*ones(N)))*VR';
end

end
