%% Import Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Data
tempData        = importdata( spec.dataName );
dataMatrix      = tempData.data;
dataHeader      = tempData.colheaders;
clear tempData

% Share Data
tempData        = importdata( spec.shareName );
shareMatrix     = tempData.data;
shareHeader     = tempData.colheaders;
clear tempData

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
    choicesetcode(index1) = bi2de(choiceset');
end

uniquecode = unique(choicesetcode);
n.choiceset = numel(uniquecode);

for k = 1:n.choiceset
    belong = choicesetcode == uniquecode(k);
    [dataR(k), data(k)] = ConstructDataGroup(dataMatrix(belong,:),n,spec);
end
