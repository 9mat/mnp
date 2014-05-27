%% Import Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Data
tempData        = importdata( spec.dataName );
dataMatrix      = tempData.data;
dataHeader      = tempData.colheaders;
clear tempData

% Share data
tempData        = importdata( spec.shareName );
shareMatrix     = tempData.data;
clear tempData;

% remove data with only 1 alternative
dataMatrix(dataMatrix(:,3) < 2, :) = [];

% for the 2160 data, station 141, 172 noone chooses the first alternative,
% which is the base alternative, thus here we simply throw them away
% dataMatrix(dataMatrix(:,1) == 141, :) = [];
% dataMatrix(dataMatrix(:,1) == 172, :) = [];
% shareMatrix(shareMatrix(:,1) == 141, :) = [];
% shareMatrix(shareMatrix(:,1) == 172, :) = [];

% ignore midgrade ethanol
if ~spec.include_emidgrade
    dataMatrix(dataMatrix(:,4)==4 | dataMatrix(:,5)==4,:) = [];
end

if spec.keep_treattype >=0
    tt = dataMatrix(:,end);
    dataMatrix(tt ~= spec.keep_treattype,:) = [];
end

shareHat = shareMatrix(:,2:end);
shareHat(:,spec.base) = NaN;
%shareHat(isnan(shareHat)) = [];

%% Separate data according to choice sets


conID           = dataMatrix( :, 2 );
alternative     = dataMatrix( :, 4 );
marketID        = dataMatrix( :, 1 );

[allID, indexID]  = unique(conID);
[uniqueID, sortID] = sort(allID);
indexID         = indexID(sortID);
marketIdByCon   = marketID(indexID);
n.maxChoice     = max(alternative);
n.marketCount   = histc(marketID, sort(unique(marketID)));
n.con           = numel(uniqueID);
mapConID(uniqueID) = 1:n.con;

% !!! Assumption: no redundant alternative 
choicesetcode   = zeros(size(conID));
for i = 1:n.con
    index1 = (conID == uniqueID(i));
    choiceset = false(1, n.maxChoice);
    choiceset(alternative(index1)) = true;
    choicesetcode(index1) = bin2dec(num2str(choiceset));
end

uniquecode = sort(unique(choicesetcode));
n.choiceset = numel(uniquecode);

%% Construct the dataR structure for each choice set
allmarkets = sort(unique(marketID));
n.market = numel(allmarkets);
choicesetsize = zeros(size(allmarkets));

dataS = cell(n.maxChoice, n.market);
dataR = cell(n.choiceset, 1);
data = cell(n.choiceset, 1);

for k = 1:n.choiceset
    belong = choicesetcode == uniquecode(k);
    [dataR{k}, data{k}] = ConstructDataGroup(dataMatrix(belong,:),n,spec);
    
    missing = true(1,n.maxChoice);
    missing(dec2bin(uniquecode(k)) == '1') = false;
    
%     for j = 1:n.maxChoice
%         if j ~= spec.base && ~missing(j)
%             dataMatrix(belong, 5) = j;
%             [dataS{k,j}, ~] = ConstructDataGroup(dataMatrix(belong,:),n,spec);
%             dataS{k,j}.draw = dataR{k}.draw;
%         end
%     end
end


for m = 1:n.market
    belong = marketID == allmarkets(m);
    i = find(belong, 1 , 'first');
      
    missing = true(1,n.maxChoice);
    missing(dec2bin(choicesetcode(i)) == '1') = false;
    
    for j = 1:n.maxChoice
        if j ~= spec.base && ~missing(j)
            dataMatrix(belong, 5) = j;
            [dataS{j,m}, ~] = ConstructDataGroup(dataMatrix(belong,:),n,spec);
        end
    end
end

%% Identify the set of estimable parameters for each data group

% Each data group has its own set of estimable parameters, depending on the
% choice set; e.g. if the choice set does not include ethanol, then
% beta_ethanol cannot be identified using the data in the group
%
% Thus we need to construct a vector to identify which parameters (from the
% original set of parameters) are estimable in each data group (vector pick
% below)

identifiable.delta = false(n.maxChoice, n.market);
for k = 1:numel(allmarkets)
    allchoices = unique(alternative(marketID == allmarkets(k)));
    identifiable.delta(allchoices,k) = true;
    choicesetsize(k) = numel(allchoices);
end
identifiable.delta(spec.base,:) = 0;
deltaindex = [0 cumsum(choicesetsize-1)'] + 1;

% Mask the identifiable parameters of the whole model (e.g. mark
% base-alternative $\beta$'s as not identifiable)
identifiable.beta_1 = true(n.maxChoice, n.prodChar);
identifiable.beta_2 = true(n.maxChoice, n.conChar);
identifiable.S      = tril(true(n.maxChoice, n.maxChoice));

identifiable.beta_1(spec.base,:)    = false;
identifiable.beta_2(spec.base,:)    = false;
identifiable.S(spec.base,:)         = false;
identifiable.S(:,spec.base)         = false;
identifiable.S(spec.scale, spec.scale) = false;

% Count the number of estimable parameters of the full model
params = fieldnames(identifiable);
for i = 1: numel(params)
    n.(params{i}) = sum(identifiable.(params{i})(:));
    paramsid.(params{i}) = zeros(size(identifiable.(params{i})));
    paramsid.(params{i})(identifiable.(params{i})) = 1:n.(params{i});
end

n.beta      = n.beta_1 + n.beta_2 + 1 + n.conGroup;
n.theta     = n.beta + n.S + n.delta;
n.maxChoice = max(dataMatrix(:,4));

paramsid.beta_1 = paramsid.beta_1 + 1 + n.conGroup;
paramsid.beta_2 = paramsid.beta_2 + 1 + n.conGroup + n.beta_1;
paramsid.delta  = paramsid.delta  + n.beta;
paramsid.S      = paramsid.S      + n.beta + n.delta;


% Index the estimable parameters in a running order

for k = 1:n.choiceset
    % Decode the choice set code to know which alternatives are missing
    missing = true(1,n.maxChoice);
    missing(dec2bin(uniquecode(k)) == '1') = false;
    
    % begin with all identifiable params of the full model
    k_identifiable = identifiable;
    
    % mark the parameters of the missing choices as unidentified
    k_identifiable.beta_1(missing,:)   = false;
    k_identifiable.beta_2(missing,:)   = false;
    k_identifiable.S(missing,:)        = false;
    k_identifiable.S(:,missing)        = false;
    
    % mark the fixed effects of outside markets as unidentified
    submarket = unique(marketID(choicesetcode == uniquecode(k)));
    dv_submarket = ismember(allmarkets, submarket);
    k_identifiable.delta = false(size(identifiable.delta));
    k_identifiable.delta(:, dv_submarket) = identifiable.delta(:,dv_submarket);

    % combine all the identified params together
    pick = sort([ (1:1+n.conGroup)'; ...
        paramsid.beta_1(k_identifiable.beta_1); ...
        paramsid.beta_2(k_identifiable.beta_2); ...
        paramsid.delta(k_identifiable.delta); ...
        paramsid.S(k_identifiable.S) ]);
        
    dataR{k}.identifiable = k_identifiable;
    dataR{k}.pick = pick;
    dataR{k}.pick_delta = sort(paramsid.delta(k_identifiable.delta));
    
    for i = 1:numel(submarket)
        m = allmarkets == submarket(i);
        k_identifiable.delta = false(size(identifiable.delta));
        k_identifiable.delta(:,m) = identifiable.delta(:,m);
        
        pick = sort([ (1:1+n.conGroup); ...
            paramsid.beta_1(k_identifiable.beta_1); ...
            paramsid.beta_2(k_identifiable.beta_2); ...
            paramsid.delta(k_identifiable.delta); ...
            paramsid.S(k_identifiable.S) ]);
        
        for j = 1:n.maxChoice
            if ~isempty(dataS{j,m})
                dataS{j,m}.identifiable = k_identifiable;
                dataS{j,m}.pick = pick;
                dataS{j,m}.pick_delta = sort(paramsid.delta(k_identifiable.delta));
            end
        end
    end
end

%% Calculate means to be used for marginal effects
meanData = ones(n.maxChoice, numel(spec.paramType));

type = mod(spec.paramType, 10);
for i = 1:size(meanData,2)
    if type(i) == 1 || type(i) == 2
        meanData(:,i) = mean(dataMatrix(dataMatrix(:,4)==spec.base,i));
    elseif type(i) == 3 || type(i) == 4
        for j = 1:n.maxChoice
            meanData(j,i) = mean(dataMatrix(dataMatrix(:,4) == j,i));
        end
    end
    
    meanData(:,3) = n.maxChoice; % choice set size
    meanData(:,4) = 1:n.maxChoice; % alternative
end
n.mfx = sum(type > 0) + (n.maxChoice-1)*sum(type > 2);


