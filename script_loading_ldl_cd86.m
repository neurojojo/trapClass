myfiles = dir('C:\#BackupMicData\**\*.mat');
myfiles = struct2table(myfiles);

[~,b] = findgroups( myfiles.folder );

myfolders = b(contains(b,'Ch2'));
%% Lining up FRET traces with UTrack traces
mytbl.tracksTable.traceLen = arrayfun(@(x) numel(x{1}), mytbl.tracksTable.x );
mytbl.tracksTable.Properties.VariableNames{2} = 'startOfTrace';

%% From fretTraces.Ch2.traceMetadata
meta_table=struct2table(fretTraces.Ch2.traceMetadata);
[x0,y0] = deal([],[]);
for i = 1:size(meta_table,1)
    x0(i) = fretTraces.Ch2.x(i,meta_table(i,:).startOfTrace);
    y0(i) = fretTraces.Ch2.y(i,meta_table(i,:).startOfTrace);
end

meta_table.x0 = x0';
meta_table.y0 = y0';
%x0 = arrayfun(@(x) x{1}(1), mytbl.tracksTable.x );

tbl1 = table( meta_table.x0,meta_table.y0 )
tbl2 = table( arrayfun(@(x) x{1}(1), mytbl.tracksTable.x ), arrayfun(@(x) x{1}(1), mytbl.tracksTable.y ) );

tbl2.id = [1:size(tbl2,1)]'

%%

consmean = @(x,min_) mean(x(x>min_));

x_ = arrayfun( @(i) consmean( fretTraces.Ch2.x(i,:),1 ), [1:size(fretTraces.Ch2.x,1)] );
y_ = arrayfun( @(i) consmean( fretTraces.Ch2.y(i,:),1 ), [1:size(fretTraces.Ch2.y,1)] );

t_ = arrayfun( @(i) fretTraces.Ch2.traceMetadata(i).startOfTrace, [1:size(fretTraces.Ch2.x,1)] );
%%

x__ = arrayfun( @(x) mean(mytbl.tracksTable(x,:).x{1}), [1:size(mytbl.tracksTable,1)] );
y__ = arrayfun( @(x) mean(mytbl.tracksTable(x,:).y{1}), [1:size(mytbl.tracksTable,1)] );
start_ = arrayfun( @(x) mean(mytbl.tracksTable(x,:).trackStart), [1:size(mytbl.tracksTable,1)] );


%%

