%%%%create graph within school
tbl = readtable('37070-0001-Data.csv');
var = tbl.Properties.VariableNames';

SCHID = table2array(tbl(:,1));
UsID = unique(SCHID);
UID = table2array(tbl(:,2));
ID = str2double(table2array(tbl(:,5)));

%spend time with: ST1-ST10
ST = table2array(tbl(:,66:75));
ST(ST==-99|ST==-97|ST==-96|ST==-95|ST==-77|ST==-55) = nan;

%change ID within school to UID
%999 as an edge recode to a newly created unique ID 
%999 as an node using the corresponding UID
ST_recode = cell(length(UsID),2);
for i = 1:length(UsID)%school ID
    t_id = ID(SCHID==UsID(i));
    t_uid = UID(SCHID==UsID(i));
    t_ST = ST(SCHID==UsID(i),:);
    for j = 1:length(SCHID(SCHID==UsID(i)))
        %first replace 999
        ind_999 = find(t_ST(j,:)==999);
        if isempty(ind_999)~=1
            t_ST(j,ind_999) = ind_999+UsID(i)*110000+j+999;
        end
        [lia,loc] = ismember(t_ST(j,:),t_id);
        t_ST(j,lia) = t_uid(loc(lia));
        %find edge not in UID & ~=999
        ind_out = find(t_ST(j,:)<999);
        if isempty(ind_out)~=1
            t_ST(j,ind_out) = UsID(i)*100000+t_ST(j,ind_out);
        end
    end
    ST_recode{i,1} = t_uid;
    ST_recode{i,2} = t_ST;
end

SimpGraph_SCH = cell(length(UsID),1);
DirGraph_SCH = cell(length(UsID),1);
for i = 1:length(UsID)
    t_uid = ST_recode{i,1};
    t_ST = ST_recode{i,2};
    edgeL = [];
    for j = 1:size(ST,2)
        edgeP = [t_uid,t_ST(:,j)];
        edgeL = [edgeL;edgeP];
    end
    %Node
    node = [edgeL(:,1);edgeL(:,2)];
    node(isnan(node)) = [];
    node = cellstr(num2str(unique(node)));
    node = strrep(node,' ','');
    %edge
    source = cellstr(num2str(edgeL(sum(isnan(edgeL),2)==0,1)));
    source = strrep(source,' ','');
    target = cellstr(num2str(edgeL(sum(isnan(edgeL),2)==0,2)));
    target = strrep(target,' ','');
    %table
    NodeTable = table(node,'VariableNames',{'Name'});
    EdgeTable = table([source,target],'VariableNames',{'EndNodes'});
    %Graph
    G_und = graph(EdgeTable,NodeTable);
    SimpGraph_SCH{i} = simplify(G_und);
    G_dir = digraph(EdgeTable,NodeTable);
    DirGraph_SCH{i} = simplify(G_dir);
end

%check node across schools is unique
Nod = [];
for i = 1:length(UsID)
    nd = str2double(table2array(SimpGraph_SCH{i}.Nodes));
    Nod = [Nod;nd];
end
%length(unique(Nod))==length(Nod)

save st_graphData_raw ST_recode SimpGraph_SCH DirGraph_SCH;
%%
%undirected degree,indegree,outdegree
Deg_und_st = cell(size(SimpGraph_SCH,1),2);
Deg_dir_st = cell(size(SimpGraph_SCH,1),3);
for i = 1:size(SimpGraph_SCH,1)
    G_und = SimpGraph_SCH{i};
    G_dir = DirGraph_SCH{i};
    %degree in undirected graph
    Deg_und_st{i,1} = str2double(table2array(G_und.Nodes));
    Deg_und_st{i,2} = centrality(G_und,'degree');
    %Indegree
    Deg_dir_st{i,1} = str2double(table2array(G_dir.Nodes));
    Deg_dir_st{i,2} = centrality(G_dir,'indegree');
    %Outdegree
    Deg_dir_st{i,3} = centrality(G_dir,'outdegree');
end

%combine
Deg1 = [];
Deg2 = [];
for i = 1:size(SimpGraph_SCH,1)
    Uid1 = Deg_und_st{i,1};
    Uid2 = Deg_dir_st{i,1};
    deg1 = Deg_und_st{i,2};
    indeg = Deg_dir_st{i,2};
    outdeg = Deg_dir_st{i,3};
    deg_1 = [Uid1,deg1];
    deg_2 = [Uid2,indeg,outdeg];
    Deg1 = [Deg1;deg_1];
    Deg2 = [Deg2;deg_2];
end
isequal(Deg1(:,1),Deg2(:,1))
Degree_st = [Deg1,Deg2(:,2:end)];

%reciporal degree
UID = [];
ST = [];
for i = 1:size(ST_recode,1)
    UID = [UID;ST_recode{i,1}];
    ST = [ST;ST_recode{i,2}];
end
edgeL = [];
for i = 1:size(ST,2)
    edgeP = [UID,ST(:,i)];
    edgeL = [edgeL;edgeP];
end
%conducted undirected graph
node = [edgeL(:,1);edgeL(:,2)];
node(isnan(node)) = [];
node = unique(node);
node = cellstr(num2str(node));%25015
node = strrep(node,' ','');

edgeL(sum(isnan(edgeL),2)~=0,:) = [];%183821
source = cellstr(num2str(edgeL(:,1)));
source = strrep(source,' ','');
target = cellstr(num2str(edgeL(:,2)));
target = strrep(target,' ','');
NodeTable = table(node,'VariableNames',{'Name'});
EdgeTable = table([source,target],'VariableNames',{'EndNodes'});

G_und = graph(EdgeTable,NodeTable);
%plot(G_und)

Edge = G_und.Edges;
Edge = str2double(table2array(Edge));
Node = G_und.Nodes;
Node = str2double(table2array(Node));
%[B,ia,ib] = unique(Edge,'rows');

deg_recip = zeros(length(Node),1);
for i = 1:length(Node)
    tp1 = find(Edge(:,1)==Node(i));
    tp2 = find(Edge(:,2)==Node(i));
    if isempty(tp1) && isempty(tp2)
        deg_recip(i) = 0;
    else
        e1 = Edge(tp1,:);
        e2 = Edge(tp2,:);
        [~,ia1,ib1] = unique(e1,'rows');
        n1 = length(ib1)-length(ia1);
        [~,ia2,ib2] = unique(e2,'rows');
        n2 = length(ib2)-length(ia2);
        deg_recip(i) = n1+n2;
    end
end

[c,ia,ib] = intersect(Degree_st(:,1),Node);
Degree_st_n = [c,Degree_st(ia,2:end),deg_recip(ib)];
tb_deg_st = array2table(Degree_st_n,'VariableNames',{'ID','undirected','indegree','outdegree','reciprocal'});

save DegData_st Degree_st_n tb_deg_st;
