
function [B, BetaNew]=extractHRFs(conds,numArrays,scanfiles,mask,basisfiles,savedirec,timeprior)
% This function can be used to extract the haemodynamic response function
% (HRF) using both the standard FIR (Finite Impulse Response) and smooth FIR. 
%
% To use this function you must input the information below:
%
% 1) conds = the conditions of interest; ie: the cond.mat files defined in the EEG-fMRI
% 2) numArrays = the number of sessions
% 3) scanfiles = the preprocessed files
% 4) mask = defined mask or region of interest
% 5) basisfiles = confound specification
% 6) savedirec = the directory which you would like to save your output
% 
% The output you get is 'Standard_FIR_Betas' with variable B (which is using 
% the standard FIR) and 'Smooth_FIR_Betas' with variable BetaNew (using the 
% smooth FIR). The size of both of these variables should be the number
% of regressors in your design matrix by the number of voxels in your mask
% (regressors x voxels).
% 
% So if you want to plot your results keep the design matrix in mind... 
% 
% An example is displayed below:
%
% conds = ['C:\data\patients\Sub_01\MRI\Analysis\session_1.cond.mat';...
%     'C:\data\patients\Sub_01\MRI\Analysis\session_2.cond.mat';...
%     'C:\data\patients\Sub_01\MRI\Analysis\session_3.cond.mat';...
%     'C:\data\patients\Sub_01\MRI\Analysis\session_4.cond.mat'];
% 
% numArrays = 4;
% 
% scanfiles = ['C:\data\patients\Sub_01\MRI\session_1\rfmri';...       
%              'C:\data\patients\Sub_01\MRI\session_2\rfmri';
%              'C:\data\patients\Sub_01\MRI\session_3\rfmri';
%              'C:\data\patients\Sub_01\MRI\session_4\rfmri'];
%          
% mask = 'C:\data\patients\Sub_01\MRI\Analysis_FIACH\mask.img';
% 
% basisfiles = ['C:\data\patients\Sub_01\MRI\session_1\rfmri\noise_basis6.txt';...
%          'C:\data\patients\Sub_01\MRI\session_2\rfmri\noise_basis6.txt';...
%          'C:\data\patients\Sub_01\MRI\session_3\rfmri\noise_basis6.txt';...
%          'C:\data\patients\Sub_01\MRI\session_4\rfmri\noise_basis6.txt'];
%      
% savedirec = 'C:\data\FIR\Patients\Sub_01';
% 
% timeprior = 0;
%
% extractHRFs(conds,numArrays,scanfiles,mask,basisfiles,savedirec,timeprior)


set(0,'DefaultFigureWindowStyle','docked')
conds = conds;
%%
freq = 128;
nsecs = 24;
%%
% How to get your data (Y)
numArrays = numArrays; % number of sessions
direc = cell(numArrays,1);
basis = cell(numArrays,1);
X_test = cell(numArrays,1);
X = cell(numArrays,1);
Y = cell(numArrays,1);
timepts = cell(numArrays,1);
new_size = cell(numArrays,1);

try
    mask = load_nii(mask);
catch
    mask = load_untouch_nii(mask);
    disp('*** load_nii did not work, using load_untouch_nii instead ***')
end
        
mask = mask.img;
onsetss = cell(numArrays,1);

for kk = 1:numArrays
    load(conds(kk,:));
    onsets = cellfun(@(x) x*2.16,onsets,'un',0);
    tr = 2.16;
    direc = scanfiles;
    files = spm_select('FPlist',direc(kk,:),'^safilt_r.*\.nii$'); % can't use the normalised files
    hdrInfo = spm_vol(files);
    data = spm_read_vols(hdrInfo);
    timepts{kk} = size(data,4);
    
    res = NaN(size(data,4), size(data,1)*size(data,2)*size(data,3));
    for iK = 1:size(data,4)
        ia = data(:,:,:,iK);
        a3 = ia(:)';
        res(iK, :) = a3;
    end
    blop = res;

    maskvec = reshape(mask,[1 size(data,1)*size(data,2)*size(data,3)]);
    maskFinal = maskvec/max(maskvec);
    Ynew = blop(:,maskFinal==1);
    Y{kk} = Ynew;
    
    %%
    % Confound specification
    
    basisfiles = basisfiles;
    basis{kk} = dlmread(basisfiles(kk,:));
    
    %%
    % Effects of Interest
    hrfLen = round(nsecs/tr)+timeprior; % FIX hrfLen-1 to account for onsets <1
    nCond = size(onsets,2);      % assuming that the first spike type is the most interesting
    onsetss{kk} = onsets;
    
    new_size{kk} = abs(min(cell2mat(onsets))-timeprior)+1;
%     new_size{kk} = abs(onsets{1}(1)-timeprior)+1; % NEEDS CHECKING
    X_FIR = zeros(timepts{kk}+round(new_size{kk}),hrfLen*nCond);
    for iC = 1:nCond
        window_onsets = onsets{iC}(:);
        onsets_FIR = window_onsets-timeprior;
        idxCols = (iC-1)*hrfLen+1:iC*hrfLen; % ie 1:21 or 1:11
        for jO = 1:numel(onsets_FIR) % 21 - number of onsets
            idxRows = onsets_FIR(jO):(onsets_FIR(jO)+hrfLen-1); %why-1?
            for kR = 1:numel(idxRows); %hrfLen
                %             if idxRows(kR)+new_size<1 %% for when you tried to index 0
                % X_FIR(1,idxCols(kR)) = 1;
                %                 continue
                %             else
                X_FIR(round(idxRows(kR)+new_size{kk}),idxCols(kR)) = 1; %278,22
            end
        end
    end

X_test{kk} = X_FIR(1:size(X_FIR,1),:);
end

% Cut from the top
for bla = 1:numArrays
yup = size(X_test{bla},1);
X{bla} = X_test{bla}(yup-timepts{bla}+1:end,:); % cutting the design matrix
end

nSess = size(X,1);

nScans = 0;
for i = 1:nSess
    nScans =  nScans + size(X{i},1);
end

nEOI = 0;
for i = 1:nSess 
    nEOI =  nEOI + size(X{i},2);
end

nRegs = nEOI + 12*nSess;

desmat = zeros(nScans,nRegs);
ndesI = 0;
pdesI = 0;
meanss = zeros(nScans,nSess);
for i =1:nSess
    desI = [X{i} basis{i}];
    meanss((ndesI+1):(ndesI+size(X{i,1})),i) = 1;
    desmat((ndesI+1):(ndesI+size(X{i,1})),(pdesI+1):(pdesI+size(desI,2))) = desI;
    ndesI = size(desI,1)+ndesI;
    pdesI = size(desI,2)+pdesI;
end
 
desmat = [desmat meanss];

Xdisp = desmat;
sdX = std(desmat); 
meanX = mean(desmat); 

for ii=1:numArrays
  Xdisp(ii,:)=(Xdisp(ii,:)-meanX);
end

for ii=1:size(desmat,2)
  Xdisp(:,ii)=(Xdisp(:,ii)./sdX(1,ii));
end

testing = find(isnan(Xdisp));
if isempty(testing)==0
    Xdisp(:,~any(~isnan(Xdisp),1)) = [];
end

figure(1)
imagesc(Xdisp)
colormap(gray)
title('Design Matrix')

%%
% Standard FIR Modeling
Y = cell2mat(Y);
B = pinv(desmat)*Y; % aka inv(X'X)*X'Y
cd(savedirec)
save('Standard_FIR_Betas_10s_test','B')

%%
% BUILDING BLOCK FOR COVARIANCE

% h = 0.09257061; % try multiplying by 1 or 2 etc.
h = 0.1; %0.001 %0.5;
v = 1;

A = 0:hrfLen+1; %12
j = repmat(A,hrfLen+2,1); %13
x = (j-j').^2;
d = hrfLen+2; %13
I = eye(d)*20; % additive white noise of variance sigma sq = 400, therefore this is 20 (Goutte et al., 2000)
E = v.*exp((-h/2).*(x)); % +I;

%% ACTUAL COVARIANCE
Sigma = zeros(size(desmat'*desmat));
seq = 1:hrfLen;
for i = 1:nSess
    R = inv(E);
    R = R(2:hrfLen+1,2:hrfLen+1);
    Sigma(seq+(i-1)*hrfLen,seq+(i-1)*hrfLen) = R;
end

% so that it doesn't matter the number of spike types you have in any sess
if nSess == 1
Sigma = (blkdiag(kron(eye(size(onsetss{1},2)),R),zeros(12),zeros(nSess))); % need to make for variables number of spike types
end

if nSess == 2
Sigma = (blkdiag(kron(eye(size(onsetss{1},2)),R),zeros(12),kron(eye(size(onsetss{2},2)),R),zeros(12),zeros(nSess))); % need to make for variables number of spike types
end

if nSess == 3
Sigma = (blkdiag(kron(eye(size(onsetss{1},2)),R),zeros(12),kron(eye(size(onsetss{2},2)),R),zeros(12),kron(eye(size(onsetss{3},2)),R),zeros(12),zeros(nSess))); % need to make for variables number of spike types
end

if nSess == 4
Sigma = (blkdiag(kron(eye(size(onsetss{1},2)),R),zeros(12),kron(eye(size(onsetss{2},2)),R),zeros(12),kron(eye(size(onsetss{3},2)),R),zeros(12),kron(eye(size(onsetss{4},2)),R),zeros(12),zeros(nSess))); % need to make for variables number of spike types
end

%% SMOOTH FIR
BetaNew = inv(desmat'*desmat + Sigma)*desmat'*Y;

cd(savedirec)
save('Smooth_FIR_Betas_10s_test','BetaNew')

% Plotting your results
figure(2)
hax=axes; 
hold on 
plot(B(1:hrfLen,:)) 
plot(median(B(1:hrfLen,:),2),'LineWidth',7)
SP=timeprior; %your point goes here 
line([SP SP],get(hax,'YLim'),'LineStyle','--','Color',[0 0 0],'LineWidth',5)
title('Standard HRF: First Spike Type')
xlim([1 hrfLen])
hold off

figure(3)
hax=axes; 
hold on 
plot(BetaNew(1:hrfLen,:)) 
plot(median(BetaNew(1:hrfLen,:),2),'LineWidth',7)
SP=timeprior; %your point goes here 
line([SP SP],get(hax,'YLim'),'LineStyle','--','Color',[0 0 0],'LineWidth',5)
title('Smooth HRF: First Spike Type')
xlim([1 hrfLen])
hold off

figure(4)
meanvoxels_standard = mean(B(1:hrfLen,:),2); 
plot(meanvoxels_standard)
SP=timeprior; %your point goes here 
line([SP SP],get(hax,'YLim'),'LineStyle','--','Color',[0 0 0],'LineWidth',3)
title('Standard HRF: First Spike Type (Over all voxels)')
xlim([1 hrfLen])

figure(5)
meanvoxels_smooth = mean(BetaNew(1:hrfLen,:),2); 
plot(meanvoxels_smooth)
SP=timeprior; %your point goes here 
line([SP SP],get(hax,'YLim'),'LineStyle','--','Color',[0 0 0],'LineWidth',3)
title('Smooth HRF: First Spike Type (Over all voxels)')
xlim([1 hrfLen])

figure(6) % not converted properly into seconds
plot([1:length(meanvoxels_smooth)]*tr-(timeprior*2.16), meanvoxels_smooth);
SP=0; %your point goes here 
line([SP SP],get(hax,'YLim'),'LineStyle','--','Color',[0 0 0],'LineWidth',3)
title('Smooth HRF: First Spike Type (Over all voxels)')
xlabel('Time (seconds)')
end
