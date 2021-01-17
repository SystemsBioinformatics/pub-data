function [required, alternatives, all_sets] = EMAF(model, constraints, exchangeRxns)

%% first transform your model to irreversible if needed
posExchange = getPosOfElementsInArray(exchangeRxns,model.rxns);
if any(model.lb(posExchange)<0)
    [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model,'sRxns',exchangeRxns,'OrderReactions',true);
    doubleDirection = cellfun(@length, rev2irrev(posExchange))==2;
    posExchangeRxns = posExchange;
    posExchangeRxns(find(doubleDirection==0)) = cell2mat(rev2irrev(posExchange(find(doubleDirection==0))));
    posExchangeRxns(find(doubleDirection==1)) = getPosOfElementsInArray(strcat(exchangeRxns(find(doubleDirection==1)),'_b'), modelIrrev.rxns);
    exchangeRxnsIrrev = modelIrrev.rxns(posExchangeRxns);
else
    modelIrrev = model;
    posExchangeRxns = posExchange;
    exchangeRxnsIrrev = modelIrrev.rxns(posExchangeRxns);
end


%% Write MILP problem
H = 1000;

S = modelIrrev.S;
[m,n_rxns] = size(S);
n_rxns_int = length(posExchangeRxns);
n_variables_total = n_rxns + n_rxns_int;


% Helper arrays for extracting solutions
sel_cont_sol = 1:n_rxns;
sel_int_sol = n_rxns+1:n_rxns+n_rxns_int;

% Set variable types
variableIDs = [modelIrrev.rxns; strcat('bin_',exchangeRxnsIrrev)];
vartype_bl(1:n_variables_total) = 'C';
vartype_bl(n_rxns+1:n_rxns+n_rxns_int) = 'B';

% Set bilevel objective
c_bl = zeros(n_variables_total,1);
c_bl(n_rxns+1:n_rxns+n_rxns_int,1) = ones(n_rxns_int,1);

% S*v = 0
A_bl = [S zeros(m,n_rxns_int)];
b_bl = zeros(m,1);
csense_bl(1:m) = 'E';

%lb
lb_bl = zeros(n_variables_total,1);
lb_bl(1:n_rxns) = modelIrrev.lb;
lb_bl(vartype_bl == 'B') = 0;

%ub
ub_bl = zeros(n_variables_total,1);
ub_bl(1:n_rxns) = modelIrrev.ub;
ub_bl(vartype_bl == 'B') = 1;

%constraints

pos_rxn_constraints = getPosOfElementsInArray(constraints.rxnList,modelIrrev.rxns);

A_constraints = zeros(length(constraints.rxnList),n_variables_total);
b_constraints = zeros(length(constraints.rxnList),1);
csense_constraints = repmat('E',1,length(constraints.rxnList));

for i = 1:length(pos_rxn_constraints)
    A_constraints(i,pos_rxn_constraints(i)) = 1;
    b_constraints(i) = constraints.values(i);
    csense_constraints(i) = constraints.sense(i);
end

A_bl = [A_bl;A_constraints];
b_bl = [b_bl;b_constraints];
csense_bl = [csense_bl,csense_constraints];

% v(j) <= ub(j)*y(j) j in selected
% v(j) - ub(j)*y(j) <= 0 in selected
A_int_constraints = zeros(n_rxns_int,n_variables_total);
for i = 1:n_rxns_int
    A_int_constraints(i,posExchangeRxns(i)) = H;
    A_int_constraints(i,n_rxns+i) = -H*modelIrrev.ub(posExchangeRxns(i));
end

A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int,1)];
csense_bl(end+1:end+n_rxns_int) = 'L';

% Construct problem structure
bilevelMILPproblem.A = A_bl;
bilevelMILPproblem.b = b_bl;
bilevelMILPproblem.c = c_bl;
bilevelMILPproblem.csense = csense_bl;
bilevelMILPproblem.lb = lb_bl;
bilevelMILPproblem.ub = ub_bl;
bilevelMILPproblem.vartype = vartype_bl;
bilevelMILPproblem.varIDs = variableIDs;
bilevelMILPproblem.contSolInd = sel_cont_sol;
bilevelMILPproblem.intSolInd = sel_int_sol;
bilevelMILPproblem.osense = 1;

%% solve MILP
tic
solution = solveCobraMILP(bilevelMILPproblem);


%% enumerate

solutions = {};
solutions_rxns = {};
solutions{1} = findSelectedInSolution(solution);
solutions_rxns{1} = exchangeRxnsIrrev(solutions{1});
lengthSolutions = length(solutions{1});

A_bl = [A_bl; [zeros(1,n_rxns), ones(1,n_rxns_int)]];
b_bl = [b_bl;lengthSolutions];
csense_bl = [csense_bl,'E'];

while strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
    pos_bin_variables = sel_int_sol(solutions{end});
    
    A_newSolution = zeros(1,n_rxns+n_rxns_int);
    A_newSolution(pos_bin_variables) = ones(1,length(pos_bin_variables));
    
    A_bl = [A_bl; A_newSolution];
    b_bl = [b_bl;lengthSolutions-1];
    csense_bl = [csense_bl,'L'];
    
    bilevelMILPproblem.A = A_bl;
    bilevelMILPproblem.b = b_bl;
    bilevelMILPproblem.csense = csense_bl;
    solution = solveCobraMILP(bilevelMILPproblem);
    if strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
        solutions{end+1} = findSelectedInSolution(solution);
        save('solutions_EMAF','solutions')
        
        solutions_rxns{end+1} = exchangeRxnsIrrev(solutions{end});
        save('solutions_rxns_EMAF','solutions_rxns')
        
%         disp('new solution found')
%         disp(length(solutions))
    end

end
toc
disp('number of solutions found:')
disp(length(solutions))
%% find required and alternative
required = solutions_rxns{1};
for i = 2:length(solutions_rxns)
    required = intersect(required, solutions_rxns{i});
end

rests = cell(length(solutions_rxns),1);
for i = 1:length(solutions_rxns)
    rests{i} = setdiff(solutions_rxns{i}, required);
end

n_alternative_groups = length(solutions_rxns{1})-length(required);
alternatives = cell(n_alternative_groups,1);

n = 1;
for i = 1:length(rests{n})
    rest_i = rests{n}(i);
    others = rests{n}(setdiff(1:length(rests{n}),i));
    alternatives_i = cellfun(@(y) setdiff(y,others), rests(find(cellfun(@(x) all(ismember(others, x)), rests))));
    alternatives{i} = union(rest_i, alternatives_i);
end
alternatives = cellfun(@(x) strjoin(x,', '), alternatives, 'UniformOutput', false);

all_sets = solutions_rxns;
end

function pos_selected = findSelectedInSolution(solution)

tol = 10^-1;
pos_selected = find(solution.int>1-tol);

end