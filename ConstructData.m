%% Import Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Data
tempData        = importdata( spec.dataName );
dataMatrix      = tempData.data;
dataHeader      = tempData.colheaders;
clear tempData

% remove data with only 1 alternative
dataMatrix(dataMatrix(:,3) < 2, :) = [];

%% Separate data according to choice sets

conID           = dataMatrix( :, 2 );
alternative     = dataMatrix( :, 4 );

uniqueID        = sort(unique(conID));
n.maxChoice     = max(alternative);
n.con           = numel(uniqueID);

choicesetcode   = zeros(size(conID));

for i = 1:n.con
    index1 = (conID == uniqueID(i));
    choiceset = false(1, n.maxChoice);
    choiceset(alternative(index1)) = true;
    choicesetcode(index1) = bi2de(choiceset);
end

uniquecode = unique(choicesetcode);
n.choiceset = numel(uniquecode);

for k = 1:n.choiceset
    belong = choicesetcode == uniquecode(k);
    [dataR(k), data(k)] = ConstructDataGroup(dataMatrix(belong,:),n,spec);
end

mask.beta_1 = ones(n.maxChoice, n.prodChar);
mask.beta_2 = ones(n.maxChoice, n.conChar);
mask.S = tril(ones(n.maxChoice, n.maxChoice));

mask.beta_2(spec.base,:) = 0;
mask.beta_1(spec.base,:) = 0;
mask.S(spec.base,:) = 0;
mask.S(:,spec.base) = 0;
spec.scale = 1 + (spec.base == 1); % if base = 1, scale = 2
mask.S(spec.scale, spec.scale) = 0;

mask.beta_1(mask.beta_1 == 1) = 1:sum(mask.beta_1(:));
mask.beta_2(mask.beta_2(:) == 1) = 1:sum(mask.beta_2(:));
mask.S(mask.S == 1) = 1:sum(mask.S(:));

for k = 1:n.choiceset
    missing = ~de2bi(uniquecode(k));
    
    pick.beta_1 = mask.beta_1;
    pick.beta_2 = mask.beta_2;
    pick.S = mask.S;
    
    pick.beta_1(missing,:) = 0;
    pick.beta_2(missing,:) = 0;
    pick.S(missing,:) = 0;
    pick.S(:,missing) = 0;
    
    pick.beta_1(pick.beta_1 == 0) = [];
    pick.beta_2(pick.beta_2 == 0) = [];
    pick.S(pick.S == 0) = [];
    
    pick.theta = (1:1+n.conGroup)'; % alpha_0, alpha_r
    temp = numel(pick.theta);
    
    pick.theta = [pick.theta;pick.beta_1(:)+temp];
    temp = numel(pick.theta);

    pick.theta = [pick.theta;pick.beta_2(:)+temp];
    temp = numel(pick.theta);

    pick.theta = [pick.theta;pick.S(:)+temp];
    dataR(k).pick = pick.theta;
end
