cd('~/Dropbox/Multitask Thesis');
addpath('~/Dropbox/Multitask Thesis/')
MMSE = csvread('ADNI1MMSE.csv',1,0);
baseline = csvread('ADNI1baseline.csv',1,0);
% in their paper, they use baseline MMSE as a predictor too!
%baseline = horzcat(baseline,MMSE(:,1));
%MMSE = MMSE(:,2:7);

train_indices = randsample(size(baseline,1),floor(size(baseline,1)*0.9));
test_indices = 1:size(baseline,1);
test_indices = test_indices(~ismember(test_indices,train_indices));
MMSEtrain = MMSE(train_indices,:);
MMSEtest = MMSE(test_indices,:);
baselinetrain = baseline(train_indices,:);
baselinetest = baseline(test_indices,:);
for i=1:5
    Xtrain{i} = baselinetrain(MMSEtrain(:,i)>=0,:);
    Ytrain{i} = MMSEtrain(MMSEtrain(:,i)>=0,i);
    Xtest{i} = baselinetest(MMSEtest(:,i)>=0,:);
    Ytest{i} = MMSEtest(MMSEtest(:,i)>=0,i);
end
% CV for parameter tuning
rho1 = exp(linspace(log(0.1),log(200),5));
rho2 = exp(linspace(log(0.1),log(200),5));
rho3 = exp(linspace(log(0.1),log(200),5));

rng(1);
cvindex = crossvalind('Kfold',size(baselinetrain,1),5);
avgerror = zeros(length(rho1),length(rho2),length(rho3));
for i=1:5 % which validation fold
    i
    Xfold_i = baselinetrain(cvindex==i,:);
    Yfold_i = MMSEtrain(cvindex==i,:);
    Xfold_others = baselinetrain(cvindex~=i,:);
    Yfold_others = MMSEtrain(cvindex~=i,:);
    for t=1:5 % extract non missing data
        X{t} = Xfold_others(Yfold_others(:,t)>=0,:);
        Xval{t} = Xfold_i(Yfold_i(:,t)>=0,:);
        Y{t} = Yfold_others(Yfold_others(:,t)>=0,t);
        Yval{t} = Yfold_i(Yfold_i(:,t)>=0,t);
    end
    for a=1:length(rho1)
        for b=1:length(rho2)
            for c=1:length(rho3)
                % disp('fitting')
                %disp(a)
                tic
                W = Least_TGL(X,Y,rho1(a),rho2(b),rho3(c));
                toc
                %                 W = Least_TGL(X,Y,rho1(a),3,4);
                for t=1:5
                    predictval{t} = Xval{t} * W(:,t);
                end
                avgerror(a,b,c) = avgerror(a,b,c) + nMSE(Yval,predictval);
            end
        end
    end
    %minparameters = minparameters + parameters;
end
% take avg, find minparam
avgerror = avgerror/5;
minerror = 100;
for a=1:length(rho1)
    for b=1:length(rho2)
        for c=1:length(rho3)
            if avgerror < minerror
                minerror = avgerror;
                minparam = [rho1(a),rho2(b),rho3(c)];
            end
        end
    end
end

W3 = Least_TGL(Xtrain,Ytrain,minparam(1),minparam(2),minparam(3));
for t=1:6
    predicted3{t} = Xtest{t} * W3(:,t);
end
nMSE(Ytest,predicted3)
wR(Ytest,predicted3)
rMSE(Ytest,predicted3)
scatter(Ytest{1},predicted3{1})
refline(1,0)
axis([10 35 10 35]);