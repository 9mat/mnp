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

%% Construct Data Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataR = ConstructDataGroup(dataMatrix,n,spec);