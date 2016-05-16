% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model deletedRxns deletedMets] = expa_elementary(model)
% 
% This function calculates a lower number of elemntary pathways for the model
% in question. This function works with the COBRA toolbox and models defined
% in the BiGG database format.
% 
% According to the paper:
% %Yeung, M., Thiele, I. and Palsson, B. O. Estimation of the number of
% %extreme pathways for metabolic networks. BMC Bioninformatics. 8:363 (2007)
% The number of extreme pathways 'explodes' as more reactions and details 
% are added to the reconstructions. The large number of such extreme
% reactions are largely due to a high number of metabolites being exchanged
% and accounted by these models, and also due to reactions which perform
% similar transformations but through different currency metabolites. This
% function simplifies the model, so to speak, in order to find more
% elementary pathways by allowing the user to remove currency metabolites
% (or any other metabolites for that matter) from the reaconstruction. Each
% elementary pathway then described has a cost, which is also returned by
% the function. The user also has the option to combine reactions that
% perform the same transformation but through different metabolites (Can go
% through nadh or nadph for example), in which case both reactions would be
% counted in the returned elementary pathway.
% 
% The function works with user defined inputs and all text output can be
% saved in a diary defined at the beginning of the function. At the end, the
% user has the option to output a video plotting all the extreme pathways
% found if a map can be given.
%
% The function also has a built in extreme pathways algorithm which
% implements the algorithm by Schilling and Palsson (2000). This function is
% rather slow, however. Part of this algorithm constitutes of building a
% matrix where all possible combinations of rows which have entries of
% opposite signs in a particular column are added in order to create a zero
% in that row. This is done for each metabolite with an unconstrained
% exchange flux. Hence, the deletion of such metabolites and their exchange
% flux as defined by the user should yield a much faster implementation of
% this function and a much lower number of extreme pathways. The outputs of
% this step of the extreme pathways calculations are also outputted in the
% command window so the user can follow the progress of this algorithm. If
% this step of the process is very slow, the user might want to consider
% deleting the metabolite hindering the process. If that is not possible,
% the user might simply take the model as it is and can use other means to
% calculate the extreme pathways. The decision of whether or not to use the
% built in expa function will determine the output of this function.
%
% Some common currency metabolites, as defined in the COBRA function
% findCarbonRxns. These are good candidates to be removed from the model:
% {'h2o','co2','o2','h2o2','nh4','no2','no3','no','h2s',...
% 'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',...
% 'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',...
% 'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf','mlthf',...
% 'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'};
% 
% INPUT:
% model - model to be modified and which extreme pathways will be calculated
% for.
% 
% OUTPUTS:
% IF USING BUILT IN EXPA FUNCTION
% model - struct array with subfield:
%   P - matrix with coefficients of extreme pathways.
%   cost - matrix of cost of each elementary parthways. Has size number of
%   pathways by length of costMets.
%   costMets - metabolites from cost. These are the metabolites that are
%   unbalanced in the matrix cost in each elementary pathway.
% IF NOT USING THE BUILT IN EXPA FUNCTION
% model - BiGG formated model, similar to inputted model, but tailored as
% the user defined. Struct fields added are:
%   del - reactions deleted since there are other, similar reactions in the
%   model. Cell array of reaction names.
%   match - cell array f reaction names of same length as del. These are
%   the reactions still in the model matching the reactions in del that
%   have been deleted.
% FOR BOTH OPTIONS
%   deletedRxns - Reactions deleted from the model by user input.
%   deletedMets - Metabolites deleted from the model according to user
%   input.
% These are both returned as a struct subfield as well if the built in expa
% function is not used.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model, deletedRxns, deletedMets] = expa_elementary(model)
%Clear objective functions to aid error
model.c = zeros(length(model.rxns),1);
fprintf('At any point, too look up a metabolite type:\nlookup {name}\n')
fprintf('At any point, too look up common currency metabolites type:\nprintall\n')
bol = askinput('Would you like to keep a diary of the command window? (y/n)',{'y' 'n'});
if bol == 'y'
    name = input('Input name of the file (no .txt extension)\n','s');
    diary([name '.txt']);
    diary on
end

for i = 1:length(model.rxns)
    if length(model.rxns{i}) > 3 && strcmp(model.rxns{i}(end-3:end),'_rev')
        model.rxns{i} = input(['Rename reaction ' model.rxns{i} ':\n'],'s');
    end
end

%Save original model
original1 = model;

%% Decompartmentalize model. This way a metabolite is completely deleted from the model.
%Remove compartment from metabolites list
for i = 1:length(model.mets)
    model.mets{i}(end-2:end) = [];
end
mets = unique(model.mets);

%Decompartmentalize model
%If a reaction is a transport reaction, however ,this step would remove
%the metabolite being trasnported from both sides and leave us with only
%currency metabolites involved in the transport. Then when we delete these
%currency metabolites, the tranport reaction would be deleted. So we
%account for that.
model.mets = cat(1,model.mets,{'transport'});
model.S = [model.S; zeros(1,length(model.rxns))];
if isfield(model,'metNames')
    model.metNames = cat(1,model.metNames,{' '});
end
if isfield(model,'metFormulas')
    model.metFormulas = cat(1,model.metFormulas,{' '});
end
if isfield(model,'metChEBIID')
    model.metChEBIID = cat(1,model.metChEBIID,{' '});
end
if isfield(model,'metKEGGID')
    model.metKEGGID = cat(1,model.metKEGGID,{' '});
end
if isfield(model,'metPubChemID')
    model.metPubChemID = cat(1,model.metPubChemID,{' '});
end
if isfield(model,'metInChIString')
    model.metInChIString = cat(1,model.metInChIString,{' '});
end
if isfield(model,'metCharge')
    model.metCharge = [model.metCharge; 0];
end
if isfield(model,'b')
    model.b = [model.b; 0];
end
for i = 1:length(mets)
    temp = find(strcmp(mets{i},model.mets));
    if length(temp)>1
       x = temp(1);
       y = temp(2:end);
       for j = length(y):-1:1
           %If we add both rows
           temp1 = model.S(x,:)+model.S(y(j),:);
           %If we subtract both rows
           temp2 = model.S(x,:)-model.S(y(j),:);
           %See if any metabolite disappears by adding the rows (appears on
           %both sides with same stoichiometric coefficient -> transport)
           temp = double(temp1~=0) - double(temp2~=0);
           model.S(x,:) = model.S(x,:)+model.S(y(j),:);
           if any(temp)
               temp = find(temp);
               for ij = 1:length(temp)
                   model.S(x,temp(ij)) = model.S(x,temp(ij))-1;
                   model.S(end,temp(ij)) = model.S(end,temp(ij))+1;
               end
           end
           model = deleteSrow(model,y(j));
       end
    end
end
for i = length(model.rxns):-1:1
    if isempty(find(model.S(:,i),1))
        model = deleteScol(model,i);
    end
end

%% Define which reactions to delete
%Removing biomass reactions is a good idea!
[t i] = max(sum(model.S~=0));
fprintf(['\n' model.rxns{i}])
bol = askinput('\nDelete any reactions initially (biomass)? (y/n)',{'y' 'n'});
if bol == 'y'
    while true
        rx = input(['Input reaction to be deleted' ...
            ' (when done type ''done'')\n'],'s');
        if strcmp(rx,'done')
            break
        elseif length(rx)>7 && isequal(rx(1:6),'lookup')
            temp = rx(8:end);
            load expa_elementary_data.mat
            pos = find(ismember(metabolites,temp));
            if ~isempty(pos)
                fprintf([formulas{pos} '    ' names{pos} '\n'])
            clear metabolites names formulas
            end
        elseif length(rx)>7 && isequal(rx,'printall')
            fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
            'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
            'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
            'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
            'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
        elseif ~ismember(rx,model.rxns)
            fprintf('WARNING: Reaction not in the model\n')
        else
            model = removeRxns(model,rx);
        end
    end
end

fin = false;
while ~fin
    while true
        fig = figure(1);
        expa_elementary_analyze(model)

        %find and print exchange reactions
        x = find(findExcRxns(model));
        fprintf(['The exchange reactions (with respective' ...
        'lower and upper bounds) are:\n'])
        for i = 1:length(x)
            fprintf([model.rxns{x(i)} '\t' num2str(model.lb(x(i))) '\t' ...
               num2str(model.ub(x(i))) '\n'])
        end

        bol = askinput(['\nWould you like to change any'...
            ' reaction bounds? (y/n)'],{'y' 'n'});
        if bol == 'y'
            while true
                rx = input('Input reaction to change bound (when done type ''done'').\n','s');
                if strcmp(rx,'done')
                    break
                elseif length(rx)>7 && isequal(rx(1:6),'lookup')
                    temp = rx(8:end);
                    load expa_elementary_data.mat
                    pos = find(ismember(metabolites,temp));
                    if ~isempty(pos)
                        fprintf([formulas{pos} '    ' names{pos} '\n'])
                    end
                    clear metabolites names formulas
                elseif length(rx)>7 && isequal(rx,'printall')
                    fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
                    'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
                    'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
                    'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
                    'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
                elseif ~ismember(rx,model.rxns)
                    fprintf('WARNING: Reaction not in the model\n')
                else
                    num = input('input new bound value:\n');
                    b = askinput(['Input bound type: upper (''u''), '...
                        'lower (''l'') or both (''b'')'],{'u' 'l' 'b'});
                    model = changeRxnBounds(model,rx,num,b);
                end
            end
        end

        %Determine if we would like to proceed
        bol = askinput(['\nWould you like to proceed with'...
            ' these metabolites and bounds? (y/n)'],{'y' 'n'});
        if bol == 'y'
            break
        end

        %Print readily removable reactions and which ones are not. Readily
        %rewmovable reactions are reactions that will not cause sinks or
        %sources if the reaction and metabolite being exchanged gets deleted
        fprintf('\nThe exchange reactions we cannot just delete are:\n')
        s = [];
        prt = 0;
        for i = 1:length(x)
            temp = full(model.S);
            y = find(temp(:,x(i)));
            %Delete metabolite and exchange reactions
            temp(:,x) = [];
            temp(y,:) = [];
            %If any reaction has either only positive or negative
            %stoichiomatric coefficients the reaction cannot be deleted
            if ~isempty(find(sum(temp > 0)==0,1)) || ...
                  ~isempty(find(sum(temp < 0)==0,1))
              fprintf([model.rxns{x(i)} '\t']);
              prt = prt+1;
                if mod(prt,5)==0 && ~(prt==0)
                    fprintf('\n')
                end
            else
                s = [s x(i)];
            end
        end
        if ~isempty(s)
            fprintf('\nThe exchange reactions we can delete are:\n')
            for i = 1:length(s)
                fprintf([model.rxns{s(i)} '\t'])
                if mod(i,5)==0
                    fprintf('\n')
                end
            end
        end
        fprintf('\n\n')
        %Delete certain metabolites or keep certain metabolites
        bol = askinput(['Would you like to delete certain exchange ' ...
            'reactions (d), keep\ncertain exchange reactions (k) or skip' ...
            ' this step (s)?'],...
            {'k' 'd' 's'});
        %If we would like to keep certain exchange reactions
        if bol=='k'
            rx = 'm';
            keep = {};
            while ~strcmp(rx,'done') && ~strcmp(rx,'return')
                rx = input(['Input exchange reaction to be kept (when ' ...
                    'done type ''done'',\nto return type ''return'')\n'],'s');
                if length(rx)>7 && isequal(rx(1:6),'lookup')
                    temp = rx(8:end);
                    load expa_elementary_data.mat
                    pos = find(ismember(metabolites,temp));
                    if ~isempty(pos)
                        fprintf([formulas{pos} '    ' names{pos} '\n'])
                    end
                    clear metabolites names formulas
                elseif length(rx)>7 && isequal(rx,'printall')
                    fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
                    'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
                    'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
                    'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
                    'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
                elseif ~ismember(rx,model.rxns) && ~strcmp(rx,'done')&& ~strcmp(rx,'return')
                    fprintf('WARNING: Reaction not in the model\n')
                elseif ~strcmp(rx,'done') && ~strcmp(rx,'return')
                    keep = cat(1,keep,rx);
                end
            end
            if ~strcmp(rx,'return')
                temp = model.rxns(findExcRxns(model));
                t = ~ismember(temp,keep);
                temp = temp(t);
                for i = 1:length(temp)
                    model = expa_elementary_deleteone(temp{i},model,original1);
                end
            end
        %If we would like to delete certain exchange reactions
        elseif bol=='d'
            while true
                rx = input(['Input exchange reaction to be deleted' ...
                    ' (when done type ''done'')\n'],'s');
                if strcmp(rx,'done')
                    break
                elseif length(rx)>7 && isequal(rx(1:6),'lookup')
                    temp = rx(8:end);
                    load expa_elementary_data.mat
                    pos = find(ismember(metabolites,temp));
                    if ~isempty(pos)
                        fprintf([formulas{pos} '    ' names{pos} '\n'])
                    end
                    clear metabolites names formulas
                elseif length(rx)>7 && isequal(rx,'printall')
                    fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
                    'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
                    'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
                    'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
                    'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
                elseif ~ismember(rx,model.rxns)
                    fprintf('WARNING: Reaction not in the model\n')
                else
                    model = expa_elementary_deleteone(rx,model,original1);
                end
            end
        end

        %Find metabolites that can be deleted. Again, these are metabolites
        %that do not cause any reactino to become a sink or a source.
        s = size(model.S,1);
        vec = {};
        for i = s:-1:1
            temp = model.S;
            temp(:,findExcRxns(model)) = [];
            temp(i,:) = [];
            if isempty(find(sum(temp>0)==0,1)) && isempty(find(sum(temp<0)==0,1))
                vec = cat(1,vec,model.mets(i));
            end
        end
        if ~isempty(vec)
            fprintf('Other metabolites that can be deleted are:\n')
            for i = 1:length(vec)
                fprintf([vec{i} '\t'])
                if mod(i,5)==0
                    fprintf('\n')
                end
            end
            fprintf('\n\n')
        end
        %Delete them if we want to
        bol = askinput(['Would you like to print all'...
            ' metabolites still in the model(list)? (y/n)'],{'y' 'n'});
        if bol == 'y'
            for i = 1:length(model.mets)
                fprintf([model.mets{i} '\t'])
                if mod(i,5)==0
                    fprintf('\n')
                end
            end
        fprintf('\n')    
        end
        bol = askinput(['Would you like to delete any additional'...
            ' metabolites? (y/n)'],{'y' 'n'});
        if bol == 'y'
            while true
                met = input(['Input metabolite to be deleted' ...
                    ' (when done type ''done'')\n'],'s');
                if strcmp(met,'done')
                    break
                elseif length(met)>7 && isequal(met(1:6),'lookup')
                    temp = met(8:end);
                    load expa_elementary_data.mat
                    pos = find(ismember(metabolites,temp));
                    if ~isempty(pos)
                        fprintf([formulas{pos} '    ' names{pos} '\n'])
                    end
                    clear metabolites names formulas
                elseif length(met)>7 && isequal(met,'printall')
                    fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
                    'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
                    'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
                    'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
                    'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
                elseif ~ismember(met,model.mets)
                    fprintf('WARNING: Metabolite not in the model\n')
                else
                    model = expa_elementary_deleteone_met(met,model,original1);
                end
            end
        end

        %Find reaction with more than one product and one reactant
        bol = askinput(['Would you like to print reactions with more ' ...
            'than one\nreactant and one product? (y/n)'],{'y' 'n'});
        if bol == 'y'
            x = ((sum(model.S>0))~=1) & ((sum(model.S>0))~=0);
            y = ((sum(model.S<0))~=1) & ((sum(model.S<0))~=0);
            t = x | y;
            printRxnFormula(model,model.rxns(t));
        end
        bol = askinput(['Would you like to delete any additional metabolites'...
            '? (y/n)'],{'y' 'n'});
        if bol == 'y'
            while true
                met = input(['Input metabolite to be deleted' ...
                    ' (when done type ''done'')\n'],'s');
                if strcmp(met,'done')
                    break
                elseif length(met)>7 && isequal(met(1:6),'lookup')
                    temp = met(8:end);
                    load expa_elementary_data.mat
                    pos = find(ismember(metabolites,temp));
                    if ~isempty(pos)
                        fprintf([formulas{pos} '    ' names{pos} '\n'])
                    end
                    clear metabolites names formulas
                elseif length(met)>7 && isequal(met,'printall')
                    fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
                    'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
                    'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
                    'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
                    'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
                elseif ~ismember(met,model.mets)
                    fprintf('WARNING: Metabolite not in the model\n')
                else
                    model = expa_elementary_deleteone_met(met,model,original1);
                end
            end
        end
        printRxnFormula(model,'ADK1');
    end
    close
    nooutput = model;
    %% Revert back to original model
    %Here we take the original model and remove the reactions and metabolites
    %we have removed above. We do this to revert compartmentalization.
    %Save deleted reactions before we delete anything
    original = original1;
    t = ~ismember(original.rxns,model.rxns);
    original.deletedRxns = original.rxns(t);

    %Remove compartment from metabolites list
    tempmets = original.mets;
    for i = 1:length(tempmets)
        tempmets{i}(end-2:end) = [];
    end
    %See which metabolites have been deleted
    t = ~ismember(tempmets,model.mets);
    %Delete them
    original.deletedMets = original.mets(t);
    %Delete them
    mets = original.mets(t);
    original = removeMetabolites(original,mets);

    %See what other reactions have been deletes
    %find transport reactions. These do not appear in the model above when we
    %removed reactions due to decompartmentalization. But we do want to keep
    %them here.
    x = ismember(original.rxns,findTransRxns(original));
    %Find reactions that have been deleted (not in the model above)
    t = ~ismember(original.rxns,model.rxns);
    %Find if reaction is not in the model and also not a transport reaction.
    t = t & ~x;
    rx = original.rxns(t);
    %Remove them.
    original = removeRxns(original,rx);
    %% Condense Reaction with the same reactants and products
    %We wish to condense reactions with same reactants and products. These
    %alternative pathways lead to a large increase in the number of extreme
    %pathways. 
    %Condense rows in S to unique rows
    [temp IA] = unique(original.S','rows');
    %find reactions that were deleted.
    deleted = find(~ismember(1:length(original.rxns),IA));
    %match reactions that were deleted with the reactions that they are equal
    %to. Here we save the name of the reaction still in the model so we can
    %assign that same name to this reactino later. Our output will then give
    %both these reactions in the extreme pathway.
    match = {};
    for i = 1:length(deleted)
        for j = 1:length(original.rxns)
            if isequal(original.S(:,deleted(i)),original.S(:,j)) && ~...
                    ismember(original.rxns{j},original.rxns(deleted))
                match = cat(2,match,original.rxns(j));
                break
            end
        end
    end
    del = original.rxns(deleted);
    %Delete deleted reactions
    original = removeRxns(original,del);

    %% Do expa
    bol = askinput(['Would you like to use the built in expa' ...
            'function? (y/n)'],{'y' 'n'});
    if bol == 'y'
        %Calculate expa
         model = expa(original,1); 
    else
        model = original;
        deletedRxns = original.deletedRxns;
        deletedMets = original.deletedMets;
        model.alternative = del;
        model.match = match;
        return
    end
    
    %Delete loops from single reversible reactions
    bol = askinput(['Would you like to proceed with expas (p) or ' ...
        'return to model tailoring (r)'],{'p' 'r'});
    if bol == 'p'
        bol = askinput('Delete loops from single reversible reactions? (y/n)',...
        {'y' 'n'});
        if bol == 'y'
            model.P = model.P(sum(double(model.P~=0),2)~=0,:);
        end
        %Restore deleted reactions to the model
        bol = askinput(['Combine reactions which go through alternative ' ...
            'currency metabolites? (y/n)'],{'y' 'n'});
        if bol == 'y'
            for i = 1:length(del)  
                %Find index
                z = findRxnIDs(model,match{i});
                %Fix it
                model.rxns = cat(1,model.rxns,del{i});
                model.P = [model.P model.P(:,strcmp(match{i},model.rxns))];
                model.lb = [model.lb; model.lb(z)];
                model.ub = [model.ub; model.ub(z)];
                model.rev = [model.rev; model.rev(z)];
                model.S = [model.S model.S(:,z)];
                model.c = [model.c; model.c(z)];
            end
        else
            for i = 1:length(del)
                %Find index of matchign reaction
                z = findRxnIDs(model,match{i});
                %Fix it
                model.P = [model.P zeros(size(model.P,1),1)];
                model.rxns = cat(1,model.rxns,del{i});
                model.lb = [model.lb; model.lb(z)];
                model.ub = [model.ub; model.ub(z)];
                model.rev = [model.rev; model.rev(z)];
                model.S = [model.S model.S(:,z)];
                model.c = [model.c; model.c(z)];
                %Find expas where matchin reaction is involved
                index = find(model.P(:,z));
                for ij = 1:length(index)
                    model.P = [model.P; model.P(index(ij),:)];
                    model.P(end,end) = model.P(end,z);
                    model.P(end,z) = 0;
                end
            end
        end
        
        [model1.rxns index] = sort(model.rxns);
        model1.P = model.P(:,index);
        model1.S = model.S(:,index);
        model1.c = model.c(index);
        model1.lb = model.lb(index);
        model1.ub = model.ub(index);
        model1.rev = model.rev(index);
        model1.mets = model.mets;
        model = model1;

        %Calculate cost of each pathway
        tempS = zeros(size(original1.S,1),size(model.S,2));
        for i = 1:length(model.rxns)
            %Find the reaction in the original model
            z = findRxnIDs(original1,model.rxns{i});
            %Copy reaction coefficients to new matrix
            tempS(:,i) = original1.S(:,z);
        end
        %Calculate cost
        cost = zeros(size(model.P,1),length(original1.mets));
        for i = 1:size(model.P,1)
            temp = tempS*model.P(i,:)';
            cost(i,:) = temp';
        end
        %Find nonempty columns
        index = find(sum(abs(cost)));
        model.cost = cost(:,index);
        model.costMets = original1.mets(index);

        deletedRxns = original.deletedRxns;
        deletedMets = original.deletedMets;

        %Arrange extreme pathways in order of number of reaction in them
        [t index] = sort(sum(model.P~=0,2));
        model.P = model.P(index,:);
        model.cost = model.cost(index,:);
        fin = true;
    else
        model = nooutput;
    end
end
%% Make video with expas
bol = askinput(['Would you like to make a video of the ' ...
        num2str(size(model.P,1)) ' extreme pathways'...
        '? (y/n)'],{'y' 'n'});
if bol == 'y'
    fprintf(['WARNING: Make sure the maps are in the same directory\n' ...
        'you are working on\n\n'])
    %Get names of files in directory
    y = dir;
    files = {};
    for i = 1:length(y)
        files = cat(1,files,y(i).name);
    end
    mapname = askinput('Input map name:',files);
    %load map
    changeCbMapOutput('matlab');
    map = readCbMap(mapname);
    bol = askinput(['Plot all reactions in the original model (a), only reactions in\n' ...
        'the ExPas (no delete reaction) (e), or all reaction in the map (m)'],{'a' 'e' 'm'});
    if bol == 'e'
        map = getrxnstp(map,model.rxns);
    elseif bol == 'a'
        map = getrxnstp(map,original1.rxns);
    end
    title = input('Write video file name (no extension)\n','s');
    hei = input('Input video heigth (pixels)\n');
    wid = input('Input video width (pixels)\n');
    fps = input('Input video speed (frames per second)\n');
    writerObj = VideoWriter([title '.avi']);
    writerObj.FrameRate = fps;
    open(writerObj);
    figure('position',[1 1 wid hei])
    model2 = model.P~=0;
    temp = full(model2(1,:))';
    drawFlux99(map,model,temp);
    frame = getframe;
    writeVideo(writerObj,frame);
    for i = 1:size(model.P,1)
        temp = full(model2(i,:))';
        drawFlux99(map,model,temp);
        frame = getframe;
        writeVideo(writerObj,frame);
    end

    close(writerObj);
end
diary off
end
%% Supporting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function asks for input using str string, and keeps asking until the
%%input given is int he cell array answers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = askinput(str,answers)
    x = input([str '\n'],'s');
    while ~ismember(x,answers)
        if length(x)>7 && isequal(x(1:6),'lookup')
            temp = x(8:end);
            load expa_elementary_data.mat
            pos = find(ismember(metabolites,temp));
            if ~isempty(pos)
                fprintf([formulas{pos} '    ' names{pos} '\n'])
            end
            clear metabolites names formulas
        elseif length(x)>7 && isequal(x,'printall')
            fprintf(['h2o, co2, o2, h2o2, nh4, no2, no3, no, h2s,\n'...
            'so3, so4, h, h2, pi, ppi, coa, accoa, ppcoa, aacoa,\n',...
            'butcoa, succoa, atp, gtp, adp, gdp, amp, gmp, nad,\n'...
            'nadp, nadh, nadph, fad, fadh, na1, ahcys, amet, thf, mlthf,\n'...
            'q8h2, q8, mql8, mqn8, 2dmmql8, 2dmmq8\n']);
        end
        x = input([str '\n'],'s');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Modified version of COBRA's drawCbMap. Does not plot text. Text Plotting
%%parts of drawCbMap have been commented out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = drawCbMap99(map,options,varargin)
%Draws a map with the specified optional parameters
%
% drawCbMap(map,options)
%
%INPUTS
% map                       COBRA map structure
%
%OPTIONAL INPUTS
% options                  Structure containing optional parameters
%   nodeWeight              Size of primary metabolite nodes  
%   nodeWeightSecondary     Size of secondary metabolite nodes
%   nodeColor               Color of metabolite nodes
%   edgeColor               Color of reaction arrows2
%   edgeArrowColor          Color of reaction arrowheads
%   edgeWeight              Width of reaction arrows
%   textSize                Font size of metabolite text
%   textColor               Text color for metaboltes
%   rxnTextSize             Font size of reaction text
%   rxnTextColor            Text color for reactions
%   otherTextColor          Color of other text
%   fileName                Name of output file
%
% varargin                  optional parameter name / parameter value pairs
%
%OUTPUT
% Displays map in a matlab figure ('matlab') or target.svg ('svg') file
% depending on value of CB_MAP_OUTPUT.
% options                  Structure containing optional parameters
%
%
% Turner Conrad     6/12/12     Added rxnTextSize default (8) to fix error 
%                               in writing reaction texts to .svg

global mapHandle
global CB_MAP_OUTPUT

%check for render output
if ~exist('CB_MAP_OUTPUT', 'var') || isempty(CB_MAP_OUTPUT)
    error('No render target specified.  Call changeCbMapOutput(outputFormat)');
end


if (nargin < 2)
    options = [];
end

if mod(length(varargin),2)==0
    for i=1:2:length(varargin)-1
        switch lower(varargin{i})
            case 'nodeweight', options.nodeWeight = cell2mat(varargin(i+1));
            case 'nodecolor', options.nodeColor = cell2mat(varargin(i+1));
            case 'edgeweight', options.edgeWeight = cell2mat(varargin(i+1));
            case 'edgecolor', options.edgeColor = cell2mat(varargin(i+1));
            case 'edgearrowcolor', options.edgeArrowColor = cell2mat(varargin(i+1));
            case 'textsize', options.textSize = cell2mat(varargin(i+1));
            case 'textcolor', options.textColor = cell2mat(varargin(i+1));
            case 'othertextsize', options.otherTextSize = cell2mat(varargin(i+1));
            case 'othertextcolor', options.otherTextColor = cell2mat(varargin(i+1));
            case 'filename', options.fileName = varargin{i+1};
            otherwise, options.(varargin{i}) = varargin{i+1};
        end
    end
else
    error('Invalid number of parameters/values');
end

%%%% Compelete the missing parts of the option
nNodes = size(map.molName,1);
nEdges = size(map.connection,1);
%Node size
if ~isfield(options,'nodeWeight')
    options.nodeWeight = ones(nNodes,1)*15;
    if strcmp(CB_MAP_OUTPUT,'svg')
        options.nodeWeight = ones(nNodes,1)*25;
    end
end

if ~isfield(options,'nodeWeightSecondary')
    options.nodeWeightSecondary = ones(nNodes,1)*10;
    if strcmp(CB_MAP_OUTPUT,'svg')
        options.nodeWeightSecondary = ones(nNodes,1)*15;
    end
end
%Node color
if ~isfield(options,'nodeColor')
    options.nodeColor = repmat([255,160,128],nNodes,1);
end
%Edge color
if ~isfield(options,'edgeColor')
    options.edgeColor = repmat([0,191,255],nEdges,1);
end
%Arrowhead color
if ~isfield(options,'edgeArrowColor')
    options.edgeArrowColor = repmat([0,0,255],nEdges,1);
end
%Edge thickness
if ~isfield(options,'edgeWeight')
    options.edgeWeight = ones(nEdges,1)*2;
    if strcmp(CB_MAP_OUTPUT,'svg')
        options.edgeWeight = ones(nEdges,1)*4;
    end
end
%Font Size
if ~isfield(options,'textSize')
    options.textSize = ones(max(nNodes,nEdges),1)*12;
    if strcmp(CB_MAP_OUTPUT,'svg')
        options.textSize = ones(max(nNodes,nEdges),1)*6;
    end
end
%Font Color
if ~isfield(options,'textColor')
    options.textColor = zeros(nNodes,3);
end

% if ~isfield(options,'otherTextSize')
%     options.otherTextSize = ones(size(map.text,1),1)*12;
%     if strcmp(CB_MAP_OUTPUT,'svg')
%         options.otherTextSize = ones(size(map.text,1),1)*135;
%     end
% end

if ~isfield(options,'otherTextColor')
    options.otherTextColor= zeros(size(map.text,1),3);
end

nodeWeight = options.nodeWeightSecondary;
nodeWeight(strcmp(map.molPrime,'Y')) = options.nodeWeight(strcmp(map.molPrime,'Y'));

if ~isfield(options,'fileName')
    options.fileName = 'target.svg';
end

if ~isfield(options,'rxnDir')
    options.rxnDir = zeros(size(map.connectionAbb,1),1);
end

if ~isfield(options,'rxnDirMultiplier')
    options.rxnDirMultiplier = 2;
end

%%%%%%%% initialization
if strcmp(CB_MAP_OUTPUT,'matlab')    % use matlab to draw the map
    clf; % this was in line 41 before
    % setting the color bar
    figure(1); colormap(cool(100))
    colorbar('location','southoutside');
    axis equal;
    hold on
elseif strcmp(CB_MAP_OUTPUT, 'java')
    % use Java/OpenGL to draw the map
    plotcbmap;
    % send the transformation coordinates
    R=map.molPosition';
    R=sort(R);
    a=max(R);
    b=min(R);
    xmax=a(1,1);
    xmin=b(1,1);
    ymax=a(1,2);
    ymin=b(1,2);
    settrans(mapHandle,xmax,xmin,ymax,ymin);
elseif strcmp(CB_MAP_OUTPUT, 'svg')
    %check fileName extension
    if isempty(regexp(lower(options.fileName),'.svg$'))
        options.fileName = strcat(options.fileName,'.svg');
    end
    textPos = (map.textPos);
    x1 = min(map.molPosition(1,:));
    if min(textPos(:,1))<x1
        x1 = min(textPos(:,1));
    end
    y1 = min(map.molPosition(2,:));
    if min(textPos(:,2))<y1
        y1 = min(textPos(:,2));
    end
    x2 = max(map.molPosition(1,:));
    if max(textPos(:,1))>x2
        x2 = max(textPos(:,1));
    end
    y2 = max(map.molPosition(2,:));
    if max(textPos(:,2))>y2
        y2 = max(textPos(:,2));
    end
    if isfield(options,'colorScale')
        numColorBins = size(options.colorScale,1);
        colorScaleWidth = 0.25*(x2-x1);
        binWidth = colorScaleWidth/numColorBins;
        colorScaleHeight = 0.05*colorScaleWidth;
        colorScaleBuffer = colorScaleHeight*1.5;
        y2 = y2+colorScaleBuffer; %add buffer for scale
    end
    if isfield(options,'fluxVarColor')
        colorScaleHeight = 0.0125*(x2-x1);
        colorScaleBuffer = colorScaleHeight*1.5;
        colorWidth = 0.025*(x2-x1);
        y2 = y2+colorScaleBuffer;
    end
    [x1,y1,x2,y2] = deal(x1-200, y1-200, x2+200, y2+200); % add buffer
    SF = .25;
    mapHandle = fopen(options.fileName, 'w');
    fprintf(mapHandle, '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
    fprintf(mapHandle,'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n');
    fprintf(mapHandle,'<svg height="%+.2f" width="%+.2f" viewBox="%+.2f %+.2f %+.2f %+.2f" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n',(y2-y1+20)*SF,(x2-x1+20)*SF,(x1-10)*SF,(y1-10)*SF,(x2-x1+0)*SF,(y2-y1+20)*SF);
    fprintf(mapHandle,'<g transform="scale(%5.3f)">\n',SF);

    %Add Scale
    if isfield(options,'overlayType')
        for i=1:numColorBins
            color = strcat('rgb(',num2str(options.colorScale(i,1)),',',num2str(options.colorScale(i,2)),',',num2str(options.colorScale(i,3)),')');
            fprintf(mapHandle,'<g style="stroke-linecap: round; stroke-linejoin: round; stroke-miterlimit: 20; fill: %s; ">\n',color);
            fprintf(mapHandle,'<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" />\n</g>\n',x2-colorScaleWidth-200*2+(binWidth*(i-1)),y2-colorScaleHeight*1.25,binWidth,colorScaleHeight);
        end
        fprintf(mapHandle, '<g style="font-family: sans; stroke: none; text-anchor:end">\n');
        fprintf(mapHandle,'<text style="fill: rgb(0,0,0); text-rendering: optimizeLegibility;" x="%.2f" y="%.2f" font-size="%dpx">%s</text>\n</g>\n',x2-colorScaleWidth-200*2-0.25*colorScaleHeight,y2-colorScaleHeight*.4,colorScaleHeight,['Scale: ' options.scaleTypeLabel '    ' options.overlayType ': ' options.overlayLB]); 
        fprintf(mapHandle, '<g style="font-family: sans; stroke: none; text-anchor:start">\n');
        fprintf(mapHandle,'<text style="fill: rgb(0,0,0); text-rendering: optimizeLegibility;" x="%.2f" y="%.2f" font-size="%dpx">%s</text>\n</g>\n',x2-200*2+0.25*colorScaleHeight,y2-colorScaleHeight*.4,colorScaleHeight,[options.overlayUB]); 
    end
    if isfield(options,'fluxVarColor')
        colorTextLabel = {'Bidirectional / reversible:', 'Unidirectional / reversible forward:', 'Unidirectional / reversible reverse:', 'Unidirectional / irreversible:'};
        color{1} = strcat('rgb(',num2str(options.fluxVarColor.biDirColor(1)),',',num2str(options.fluxVarColor.biDirColor(2)),',',num2str(options.fluxVarColor.biDirColor(3)),')');
        color{2} = strcat('rgb(',num2str(options.fluxVarColor.uniDirFwdColor(1)),',',num2str(options.fluxVarColor.uniDirFwdColor(2)),',',num2str(options.fluxVarColor.uniDirFwdColor(3)),')');
        color{3} = strcat('rgb(',num2str(options.fluxVarColor.uniDirRevColor(1)),',',num2str(options.fluxVarColor.uniDirRevColor(2)),',',num2str(options.fluxVarColor.uniDirRevColor(3)),')');
        color{4} = strcat('rgb(',num2str(options.fluxVarColor.uniDirIrrColor(1)),',',num2str(options.fluxVarColor.uniDirIrrColor(2)),',',num2str(options.fluxVarColor.uniDirIrrColor(3)),')');
        for i=2:3:11
            fprintf(mapHandle, '<g style="font-family: sans; stroke: none; text-anchor:end">\n');
            fprintf(mapHandle,'<text style="fill: rgb(0,0,0); text-rendering: optimizeLegibility;" x="%.2f" y="%.2f" font-size="%dpx">%s</text>\n</g>\n',i*(x2-x1)/12,y2-colorScaleHeight*.4,colorScaleHeight,colorTextLabel{ceil(i/3)});            
            fprintf(mapHandle,'<g style="stroke-linecap: round; stroke-linejoin: round; stroke-miterlimit: 20; fill: %s; ">\n',color{ceil(i/3)});
            fprintf(mapHandle,'<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" />\n</g>\n',i*(x2-x1)/12,y2-colorScaleHeight*1.25,colorWidth,colorScaleHeight);
        end
    end
end

%%%%% actual map drawing code
% draw other shapes
%if isfield(map,'shapeThickness')
%    for i = 1:size((map.shapeThickness),1)
%        drawShape(map.shapeType(i,1),map.shapePos(i,1:2),map.shapeSize(i,1:2),map.shapeColor(i,1:3),map.shapeThickness(i,1),map.shapeStyle(i,1));
%    end
%end
% draw the connection segments traversing through the connection matrix
for i = 1:(size((map.connection),1))
    drawLine(map.connection(i,1),map.connection(i,2),map,options.edgeColor(i,:),options.edgeArrowColor(i,:),options.edgeWeight(i),nodeWeight,options.rxnDir(i),options.rxnDirMultiplier);
end
% draw the circles representing molecules
for i = 1:size((map.molPosition),2)
    drawCircle(map.molPosition(:,i),nodeWeight(i),options.nodeColor(i,:));
end
% draw texts
%{
for i = 1:length(map.text)
    textFont =map.textFont{i};
    if regexp(textFont,'@')
        [textFont, textSize] = strtok(textFont,'@');
        textSize = str2num(regexprep(textSize,'@',''));
    elseif(map.textSize(i) >= 60)
        textSize = 60;
    else
        textSize = map.textSize(i);
    end
    if find(regexp(textFont,'Italic'))
        textStyle = 'italic;';
    else
        textStyle = '';
    end
    %textFont
    % if find(regexp(textFont,' B')) || find(regexp(textFont(end),'B'))
    if (find(regexp(textFont,'B')))
        textWeight = 'bold';
        textFont = regexprep(textFont,' B','');
    else
        textWeight = '';
    end
    if isfield(options,'otherTextSize'), textSize = options.otherTextSize(i); end
    drawText(map.textPos(i,1),map.textPos(i,2),map.text{i,1},textSize,textStyle,options.otherTextColor(i,:),lower(textFont),textWeight,true);
end
% Write Metabolite Label
for i = 1:size((map.molPosition),2)  
    % write the labels for molecules
    if(options.textSize(i) ~= 0)
        drawText(map.molLabelPos(i,1),map.molLabelPos(i,2),map.molAbbreviation{i},options.textSize(i),'',options.textColor(i,:));
    end
end
% Write Reaction Label
for i = 1:size(map.rxnLabelPosition,2)
    if ~any(isnan(map.rxnLabelPosition(:,i)))
      if isfield(options, 'rxnTextSize')
        drawText(map.rxnLabelPosition(1,i),map.rxnLabelPosition(2,i),map.connectionAbb{find(map.rxnIndex(i)==map.connection,1)},options.rxnTextSize(i),'italic');
      else
        drawText(map.rxnLabelPosition(1,i),map.rxnLabelPosition(2,i),map.connectionAbb{find(map.rxnIndex(i)==map.connection,1)},8,'italic');
      end
    end
end
%}
if strcmp(CB_MAP_OUTPUT,'matlab')
    hold off;
elseif strcmp(CB_MAP_OUTPUT,'java')

elseif strcmp(CB_MAP_OUTPUT,'svg')

    fprintf(mapHandle,'</g>\n');
    fprintf(mapHandle,'</svg>\n');
    fclose(mapHandle);
    display('Document Written')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Modified version of COBRA's drawFlux. Calls drawCbMap99.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = drawFlux99(map,model,flux,options,varargin)
%drawFlux overlays a flux distribution onto a reaction map
%
% options = drawFlux(map,model,flux,options,varargin)
%
%INPUTS
% map               map structure
% model             COBRA model structure
% flux              Flux vector to overlay
%
%OPTIONAL INPUTS
% Optional parameters can be set using either the options structure, a
% parameter name / value pair input arguments, or a combination of both.
%
% options            Structure containing optional parameters
%   lb                Lower limit to round smaller values up to.
%   ub                Upper limit to round larger values down to.
%   colorScale        Colormap
%   zeroFluxWidth     Width of arrows of reactions which carry zero flux.
%   zeroFluxColor     Color of arrows of reactions which carry zero flux.
%   fileName          Name of output file
%   rxnDirMultiplier  scaling value of arrows denoting flux direction
%
%  Note: see setMapOptions for additional options.
%
% varargin          optional parameter name / parameter value pairs
%
%OUTPUT
% options           Structure containing optional parameters.
%
%
%

if nargin<4, options=[]; end
%Parse optional parameters
if mod(length(varargin),2)==0
    for i=1:2:length(varargin)-1
        options = setMapOptions(options,map,model,varargin{i},varargin{i+1});
    end
else
    error('Invalid number of parameters/values');
end

if ~isfield(options,'colorScale')
    options.colorScale = cool(100);
end
if ~isfield(options,'scaleType'), options.scaleType=1; end
if ~isfield(options,'lb'), lb=[];else lb = options.lb; end
if ~isfield(options,'ub'), ub=[];else ub = options.ub; end
if ~isfield(options,'rxnDirMultiplier'), options.rxnDirMultiplier = 2; end
if ~isfield(options,'rxnDirFlag'), rxnDirFlag = false; else rxnDirFlag = options.rxnDirFlag; end
rxnListZero = model.rxns(abs(flux)<=1e-9);
absFlag=false;
switch lower(options.scaleType)
    case {1, 'linear'}
        options.scaleTypeLabel='Linear;';
    case {2 ,'linear absolute'}
        flux=abs(flux);
        absFlag=true;
        options.scaleTypeLabel='Linear absolute;';
    case {3,'log10'}
        flux = log10(abs(flux));
        rxnListZero = model.rxns(isinf(flux));
        options.scaleTypeLabel='Log10;';
end
if ~isempty(ub)
    flux(flux>ub)=ub;
    options.overlayUB = [num2str(ub) '+'];
    fluxMax = ub;
else
    options.overlayUB = num2str(max(flux));
    fluxMax = max(flux);
end
if ~isempty(lb)
    flux(flux<lb)=lb;
    options.overlayLB = [num2str(lb) '-'];
    fluxMin = lb;
elseif absFlag
    options.overlayLB = '0';
    fluxMin = 0;
else
    fluxMin = min(flux(~isinf(flux)));
    options.overlayLB = num2str(fluxMin);
end
if isempty(find(options.colorScale>1, 1))
    options.colorScale = round(options.colorScale*255);
end
flux2 = flux-fluxMin;
if (fluxMax-fluxMin~=0), flux2 = flux2/(fluxMax-fluxMin); end
color = getColorFromColorScale(flux2,options.colorScale);
if isfield(options,'zeroFluxWidth')
    global CB_MAP_OUTPUT
    if ~isfield(options,'edgeWeight')
        s= size(map.connection);      
        if strcmp(CB_MAP_OUTPUT,'svg')
            options.edgeWeight = ones(s(1),1)*9;
        else
            options.edgeWeight = ones(s(1),1)*2;
        end
    end
    options.edgeWeight(ismember(map.connectionAbb,rxnListZero))=options.zeroFluxWidth;
end
if isfield(options,'zeroFluxColor')
    zeroFluxRxns = find(ismember(model.rxns,rxnListZero));
    color(zeroFluxRxns,:)=repmat(options.zeroFluxColor,size(zeroFluxRxns,1),1);
end

%rxnDirectionality
if rxnDirFlag
    options.rxnDir = zeros(length(map.rxnIndex),1);
    options.rxnDir(ismember(map.connectionAbb,model.rxns(flux>0))) = 1;
    options.rxnDir(ismember(map.connectionAbb,model.rxns(flux<0))) = -1;
end



options = setMapOptions(options,map,model,'edgeColor',color);
options.colorScale=options.colorScale;
options.lb = fluxMin;
options.ub = fluxMax;
options.overlayType = 'Flux';
drawCbMap99(map,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function takes the input map and returns it only with the reactions
%%rxnstp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = getrxnstp(map,rxnstp)

%delete reactions we don't want
td = zeros(1,length(map.connectionAbb));
for i=1:length(rxnstp)
    for j = 1:length(map.connectionAbb)
        if strcmp(map.connectionAbb{j},rxnstp{i})
            td(j) = 1;
        end
    end
end
sz = 1;
while sz <= length(td)
    if td(sz)==0
        ij = find(map.connection(sz,1)==map.rxnIndex);
        map.rxnPosition(:,ij) = [];
        map.rxnLabelPosition(:,ij) = [];
        map.rxnIndex(ij) = [];
        td(sz) = [];
        map.connection(sz,:) = [];
        map.connectionAbb(sz) = [];
        map.connectionName(sz) = [];
        map.connectionReversible(sz) = [];
    else
        sz = sz+1;
    end
end
%delete other texts we don't want
map.text = [];
map.textFont = [];
map.textPos = [];
map.textSize = [];
%delete nodes we don't want
sz = size(map.connection);
clear td
td = zeros(1,length(map.molIndex));
for i = 1:sz(1)
    for j = 1:length(map.molIndex)
        if map.connection(i,1) == map.molIndex(j)
            td(j) = 1;
        end
        if map.connection(i,2) == map.molIndex(j)
            td(j) = 1;
        end
    end
end
clear sz
sz = 1;
while sz <= length(map.molIndex)
    if td(sz) == 0
        td(sz) = [];
        map.molPosition(:,sz) = [];
        map.molIndex(sz) = [];
        map.molName(sz) = [];
        map.molAbbreviation(sz) = [];
        map.molLabelPos(sz,:) = [];
        map.molPrime(sz) = [];
        map.molCompartment(sz) = [];
    else
        sz = sz+1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Analyze model for number of reactions each metabolite is on. Currency
%%metabolites should appear in more of them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expa_elementary_analyze(model)

 %Get metabolites that appear the most
    M = sum(model.S~=0,2);
    [M index] = sort(M);
    mets1 = model.mets(index);
    if length(M)>24
        subplot(2,1,2)
        bar(M(end-24:end))
        set(gca,'XTick',1:25,'XTickLabel',mets1(end-24:end),...
                        'FontSize',10)
        title('Metabolite Appearance in the Model','FontSize',16)
        xlabel('Metabolite','FontSize',14)
        ylabel('Appearances','FontSize',14)
        axis([0 26 0 1.1*max(M(end-24:end))])
    else
        subplot(2,1,2)
        bar(M)
        set(gca,'XTick',1:length(M),'XTickLabel',mets1,...
            'FontSize',10)
        title('Metabolite Appearance in the Model','FontSize',16)
        xlabel('Metabolite','FontSize',14)
        ylabel('Appearances','FontSize',14)
    end
    temp = double(model.S~=0);
    temp = full(temp);
    temp = temp*temp';
    for ij = 1:size(temp,1)
        temp(ij,ij) = 0;
    end
    [y i] = max(temp);
    [y index] = sort(y);
    mets2 = model.mets(index);
    if length(y)>24
        subplot(2,1,1)
        bar(y(end-24:end),'r')
        set(gca,'XTick',1:25,'XTickLabel',mets2(end-24:end),...
                        'FontSize',10)
        title('Metabolites Coappearance in the Model','FontSize',16)
        xlabel('Metabolite','FontSize',14)
        ylabel('Coappearances','FontSize',14)
        axis([0 26 0 1.1*max(y(end-24:end))])
    else
        subplot(2,1,1)
        bar(y,'r')
        set(gca,'XTick',1:length(y),'XTickLabel',mets2,...
            'FontSize',10)
        title('Metabolites Coappearance in the Model','FontSize',16)
        xlabel('Metabolite','FontSize',14)
        ylabel('Coappearances','FontSize',14)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Built in function to calculate the extreme pathways. This function
%%implements the algorithm in:
%%Schilling, C. H., Letscher, D. & Palsson, B. O. Theory for the systemic 
%%definition of metabolic pathways and their use in interpreting metabolic 
%%function from a pathway-oriented perspective. Journal of theoretical 
%%biology. 203, 229–48 (2000).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model fullmodel] = expa(model,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [model fullmodel] = expaths(model,print)
%
%This function finds the extreme pathways of a stoichiometric model
%according to the algorithm by Schilling, Letscher and Palsson:
%
%Schilling, C. H., Letscher, D. & Palsson, B. O. Theory for the systemic 
%definition of metabolic pathways and their use in interpreting metabolic 
%function from a pathway-oriented perspective. Journal of theoretical 
%biology 203, 229–48 (2000). 
%
%This algorithm is defined in appendix B of the citation above. Steps of
%the algorithm and relevant excerpts from the paper are quoted as comments
%throuhghout the code
%
%INPUT:
%model - model is inputed in the format used with the COBRA toolbox.
%Minimum required fields are
%   S - stoichiometric coefficient matrix
%   lb - lower bound for reactions
%   ub - upper bound for reactions
%   rev - reversibility flag for each reaction in the model. Freely
%       exchanged metabolites (From exchange reactions) are not determined to
%       be so through their reversibility flag, but through non-zero upper and
%       lower bounds. Inputs have a zero upper bound and outputs have a zero
%       lower bound
%   c - objective function
%   rxns - reaction names
%
%OPTIONAL INPUT
%Print - If this input is given and is 1, 'y' or 'Y' the function prints
%the progress of the algorithm.
%
%OUTPUTS
%fullmodel - decomposed model of original model inputed. Reversible
%   reactions get decomposed into two forward reactions. All fields mentioned
%   above get adjusted, as well as rules,rxnGeneMat,grRules,subSystems,
%   confidenceScores,rxnReferences,rxnECNumbers,rxnNotes,rxnNames. These are
%   commonly used COBRA fields. most importantly, a matrix P is added as a
%   field and corresponds to the extreme pathways.
%model - the model outputted matches the original model inputed. This is a
%   collapsed version of the 'fullmodel' output, where the reversible
%   reactions that had been decomposed are combined again. The extreme
%   pathways are combined as absolute value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 1 && (varargin{1} == 'y' || varargin{1} == 'Y' || varargin{1} == 1)
    print = 1;
else
    print = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = full(model.S);
rxns = model.rxns;
ub = model.ub;
lb = model.lb;
rev = model.rev;

%% Arrange model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Put things in order first. Rearrange reactions so that exchange fluxes are
%in the end and decompose all reversible reactions into two opposing
%fluxes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make sure vectors are row vectors
if iscolumn(rxns)
    rxns = rxns';
end
if iscolumn(lb)
    lb = lb';
end
if iscolumn(ub)
    ub = ub';
end
if iscolumn(rev)
    rev = rev';
end
%Define what kind of flux each of the fluxes are. 'x' means exchange flux
%and i means 'internal'
sz = size(S);
fluxtype = cell(1,sz(2));
xcount = 0;
icount = 0;
for i = 1:sz(2);
    if length(find(S(:,i))) == 1
        fluxtype{i} = 'x';
        xcount = xcount+1;
    else
        fluxtype{i} = 'i';
        icount = icount+1;
    end
end

%Decompose all reversible reactions into two fluxes. "All internal
%reactions that are considered to be capable of operating in a reversible
%fashion are considered as two fluxes occuring in opposite directions,
%therefore constraining all internal fluxes to be non-negative.
[m n] = size(S);
%make room
S = [S zeros(m,n)];
lb = [lb zeros(1,length(lb))];
ub = [ub zeros(1,length(ub))];
rev = [rev zeros(1,length(rev))];
fluxtype = cat(2,fluxtype,cell(1,length(fluxtype)));
rxns = cat(2,rxns,cell(1,length(rxns)));
%Decompose
count = 1;
ixs = zeros(1,2*n);
ixs(1:n) = 1:n;
for i = 1:n
    %If reaction is reversible
    if model.rev(i) && fluxtype{i} == 'i'
        %Add reverse of the reaction
        S(:,n+count) =  -S(:,i);
        fluxtype{n+count} = 'i';
        rxns{n+count} = [rxns{i} '_rev'];
        ixs(n+count) = i;
        ub(n+count) = -lb(i);
        lb(i) = 0;
        lb(n+count) = 0;
        rev(n+count) = 1;
        count = count+1;
    end
end
%Clear empty space
for i = 2*n:-1:1
    if ~isempty(find(S(:,i),1))
        i = i+1;
        S(:,i:2*n) = [];
        lb(i:2*n) = [];
        ub(i:2*n) = [];
        rev(i:2*n) = [];
        rxns(i:2*n) = [];
        fluxtype(i:2*n) = [];
        ixs(i:2*n) = [];
        break
    end
end
[~, n] = size(S);

%Arrange matrix so that exchange fluxes come after the internal fluxes
%%"We typically structure the stoichiometric matrix so that the first
%%series of columns represent the internal fluxes and the remaining columns
%%represent the exchange fluxes"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rearrange reactions by name
[rxns indexes] = sort(rxns);
ub = ub(indexes);
lb = lb(indexes);
rev = rev(indexes);
fluxtype = fluxtype(indexes);
S = S(:,indexes);
ixs = ixs(indexes);
%rearrange reactions by type
[fluxtype indexes] = sort(fluxtype);
rxns = rxns(indexes);
ub = ub(indexes);
lb = lb(indexes);
rev = rev(indexes);
S = S(:,indexes);
ixs = ixs(indexes);

%Save new model
model.S = sparse(S);
model.rxns = rxns;
model.ub = ub;
model.lb = lb;
model.fluxtype = fluxtype;
model.rev = rev;
model.c = zeros(1,n);
if isfield(model,'rules')
    model.rules = model.rules(ixs);
end
if isfield(model,'rxnGeneMat')
    model.rxnGeneMat = model.rxnGeneMat(ixs,:);
end
if isfield(model,'grRules')
    model.grRules = model.grRules(ixs);
end
if isfield(model,'subSystems')
    model.subSystems = model.subSystems(ixs);
end
if isfield(model,'confidenceScores')
    model.confidenceScores = model.confidenceScores(ixs);
end
if isfield(model,'rxnReferences')
    model.rxnReferences = model.rxnReferences(ixs);
end
if isfield(model,'rxnECNumbers')
    model.rxnECNumbers = model.rxnECNumbers(ixs);
end
if isfield(model,'rxnNotes')
    model.rxnNotes = model.rxnNotes(ixs);
end
if isfield(model,'rxnNames')
    model.rxnNames = model.rxnNames(ixs);
end

%% Begin Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Begin Extreme pathway algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = full(model.S);
sz = size(S);
%The algorithm begins with the formulation of an initial matrix consisting
%of an n by n identity matrix appended to the transpose of the
%stoichiometry matrix S, ST
[metsnum rxnnum] = size(S);
T0 = zeros(rxnnum,rxnnum + metsnum);
T0(:,1:rxnnum) = eye(rxnnum);
T0(:,rxnnum+1:end) = S';

%Then we examine the constraints on each of the exchange fluxes as given in
%the equation: alphaj <= bj <= betaj. If the exchange flux is constrained
%to be positive nothing is done. however, if the exchange flux is 
%constrained to be negative, then we multiply the corresponding row of the 
%initial matrix by -1. If the exchange flux is unconstrained then we move 
%the entire row to a temporaty matrix, TE. This completes the 
%initialization of the first tableau, TO.
remove = [];
szT0 = size(T0);
TE = [];
for i = szT0(1)-xcount:szT0(1)
    if lb(i) < 0 && ub(i) > 0
        TE = [TE; T0(i,:)];
        remove = [remove i];
    elseif lb(i) < 0 && ub(i)<=0
        T0(i,:) = T0(i,:)*(-1);
    end
end

for i = length(remove):-1:1
    T0(remove(i),:) = [];
end
szT0 = size(T0);
%Step 1: Identify all metabolites that do not have an unconstrained
%exchange flux associated with them. The total number of such metabolites
%is denoted by mu.
if ~isempty(TE)
    c = [];
    szTE = size(TE);
    test = TE(:,rxnnum+1:end);
    c = find(sum(abs(test),1)==0);
else
    c = 1:length(model.mets);
end
mu = length(c);
clear test x
T0 = sparse(T0);
TE = sparse(TE);
%Step 5: repeat teps 2 to 4 for all of the metabolites that do not have an
%unconstrained exchange flux operating on the metabolite. This is performed
%by this for loop containing steps 2 to 4.
tic
alreadysaidso = false;
while ~isempty(c)
%Determine Best c to be used min(pos*neg). By using the best argument c
%possible we minimize the computational cost.
posnum = sum(T0(:,rxnnum+c)>0);
negnum = sum(T0(:,rxnnum+c)<0);
[~ , i] = min(posnum.*negnum);
%Print progress:
if print
    fprintf([num2str(length(c)) ' metabolites left'...
        '\t' model.mets{c(i)} '\n'])
    fprintf(['Size of T is: ' num2str(size(T0,1)) '\n'])
end
tic
%Step 2: Begin forming the new matrix Tx by copying all rows from Tx-1
%which contain a zero in the column of ST that corresponds  to the first
%metabolite indentified in step 1, denoted by the index c
    %Initialize T1
    %find rows that contain a zero in the metabolite's column
    temp = find(T0(:,rxnnum+c(i))==0);
    %copy rows to T1
    T1 = T0(temp,:);
    %find all other rows
    x = 1:size(T0,1);
    y = ~ismember(x,temp);
    x = x(y);
    %clear rows from T0
    T0 = T0(x,:);
    szT0 = size(T0);
    
%Step 3: Of the remaining rows in T(x-1) add together all possible
%combinations of rows which contain values of the opposite sign in the
%colum c, such that the addition produces a zero in this column.
%Find vectors of positive and negative numbers
pos = find(T0(:,rxnnum+c(i))>0);
neg = find(T0(:,rxnnum+c(i))<0);
%Get matrix of rows with negative coefficients
T0neg = T0(neg,:);
T0neg = clearrows(T0neg);
%Get matrix of rows with positive coefficients
T0pos = T0(pos,:);
T0pos = clearrows(T0pos);
if print
fprintf([num2str(size(T0pos,1)) ' positive rows\n'])
fprintf([num2str(size(T0neg,1)) ' negative rows\n'])
end
if ~isempty(pos) && ~isempty(neg)
    if alreadysaidso
        saidsofornow = true;
    else
        saidsofornow = false;
    end
    if size(T0pos,1) < size(T0neg,1)
        %get positive coefficients
        posent = T0(pos,rxnnum+c(i));
        %Make diagonal matrix with negative coefficients
        negM = abs(diag(T0(neg,rxnnum+c(i))));
        s = size(negM,1);
        for ij = 1:length(pos)
            if saidsofornow
                fprintf(['Calculating ' num2str(ij) ' of ' num2str(length(pos)) ...
                    ', size of T is ' num2str(size(T1,1)) '\n'])
            end
            temp = repmat(T0(pos(ij),:),s,1);
            T1 = [T1; negM*temp + diag(ones(1,s)*posent(ij))*T0neg];
            %Clear conical dependence (Step 4). Auxiliar function.
            T1 = clearrows(T1);
            if toc>15 && ~saidsofornow
                fprintf(['Currentyl on ' num2str(ij) ' of ' num2str(length(pos)) ...
                    ', size of T is ' num2str(size(T1,1)) '\n'])
                fprintf('Computation is getting costly! Metabolites still to go:\n')
                for ijk = 1:length(c)
                    fprintf([model.mets{c(ijk)} '\n'])
                end
                bol = askinput('Would you like to proceed (p) or exit expa calculation (r)',{'p' 'r'});
                if bol == 'r'
                    model = []; 
                    fullmodel = [];
                    fprintf('Exiting expa calculation. Empty outputs returned\n\n')
                    return
                end
                bol = askinput(['Would you like me to keep asking whether you would like to procced (a)?\n',...
                    'or would you like to go on with expa calculation to the end (e)?'], {'a' 'e'});
                
                saidsofornow = true;
                if bol == 'e'
                    alreadysaidso = true;
                end
            end
        end
    else
        %get negative coefficients
        negent = T0(neg,rxnnum+c(i));
        %Make diagonal matrix with positive coefficients
        posM = diag(T0(pos,rxnnum+c(i)));
        s = size(posM,1);
        for ij = 1:length(neg)
            if saidsofornow
                fprintf(['Calculating ' num2str(ij) ' of ' num2str(length(neg)) ...
                    ', size of T is ' num2str(size(T1,1)) '\n'])
            end
            temp = repmat(T0(neg(ij),:),s,1);
            T1 = [T1; posM*temp + abs(diag(ones(1,s)*negent(ij)))*T0pos];
            %Clear conical dependence (Step 4). Auxiliar function.
            T1 = clearrows(T1);
            if toc>15 && ~saidsofornow
                fprintf(['Currentyl on ' num2str(ij) ' of ' num2str(length(neg)) ...
                    ', size of T is ' num2str(size(T1,1)) '\n'])
                fprintf('Computation is getting costly! Metabolites still to go:\n')
                for ijk = 1:length(c)
                    fprintf([model.mets{c(ijk)} '\n'])
                end
                bol = askinput('Would you like to proceed (p) or exit expa calculation (r)',{'p' 'r'});
                if bol == 'r'
                    model = []; 
                    fullmodel = [];
                    fprintf('Exiting expa calculation. Empty outputs returned\n\n')
                    return
                end
                bol = askinput(['Would you like me to keep asking whether you would like to procced (a)?\n',...
                    'or would you like to go on with expa calculation to the end (e)?'], {'a' 'e'});
                saidsofornow = true;
                if bol == 'e'
                    alreadysaidso = true;
                end
            end
        end
    end
end
clear T0 temp pos neg T0pos T0neg posM negM posent negent
T0 = T1;
szT0 = size(T0);
clear T1
c(i) = [];
end

T1 = T0;
clear T0
clear temp
%note that the number of extreme pathways will be equal to the number of
%rows in T1
%Step 6: Next we append TE to the bottom of T1

T1 = [T1; TE];
%Step 8: follow the same procedure as in step 7 for each of the columns on
%the right side fo the tableau containing nonzero entries

for i = sz(2)+1:size(T1,2)
    if ~isempty(find(T1(:,i),1))
%Step 7: Starting at the n+1 column (or the first non-zero column on the
%right side), if Ti,(n+1) does not equal zero, then add the corresponding
%non-zero row from TE to row i so as to produce a zero in the (n+1) column.
%This is done by simply multiplying the corresponding row in TE by Ti,(n+1)
%and adding this to row i. Repeat this procedure for each of the rows in
%the upper portion of the tableau so as to create zeros in the entire upper
%portion of the (n+1) column. When finished, remove the row in TE
%corresponding to the exchange flux for the metabolite just balanced.
        x = find(T1(:,i));
        if length(x) > 1
            for j = 1:length(x)-1
                T1(x(j),:) = T1(x(j),:) + T1(x(end),:)*T1(x(j),i);
            end
        end
        T1(x(end),:) = [];
    end
end
%The final tableau, Tfinal, wil contain the transpose of the matrix P
%containing the extreme pathways in place of the original identity matrix.
model.P = T1(:,1:sz(2));
fullmodel = model;
%Revert model back to original model
s = length(model.rxns);
for i = s:-1:1
    %If reaction is the reverse of another reaction enter the loop
    if length(model.rxns{i})>3 && strcmp(model.rxns{i}(end-3:end),'_rev')
        %find index if original reaction
        j = find(ismember(model.rxns,model.rxns{i}(1:end-4)));
        if ~all(size(model.P(:,j))==size(model.P(:,i)))
            keyboard
        end
        %Fix matrix P and upper and lower bounds
        model.P(:,j) = model.P(:,j) - model.P(:,i);
        model.lb(j) = -model.ub(i);
        %Define reaction as reversible
        model.rev(j) = 1;
        %Delete _rev reaction
        model.S(:,i) = [];
        model.rxns(i) = [];
        model.rev(i) = [];
        model.lb(i) = [];
        model.ub(i) = [];
        model.c(i) = [];
        model.P(:,i) = [];
        model.fluxtype(i) = [];
        if isfield(model,'rules')
            model.rules(i) = [];
        end
        if isfield(model,'rxnGeneMat')
            model.rxnGeneMat(i,:) = [];
        end
        if isfield(model,'grRules')
            model.grRules(i) = [];
        end
        if isfield(model,'subSystems')
            model.subSystems(i) = [];
        end
        if isfield(model,'confidenceScores')
            model.confidenceScores(i) = [];
        end
        if isfield(model,'rxnReferences')
            model.rxnReferences(i) = [];
        end
        if isfield(model,'rxnECNumbers')
            model.rxnECNumbers(i) = [];
        end
        if isfield(model,'rxnNotes')
            model.rxnNotes(i) = [];
        end
        if isfield(model,'rxnNames')
            model.rxnNames(i) = [];
        end
    end
end
    fullmodel.rev = fullmodel.rev';
    fullmodel.lb = fullmodel.lb';
    fullmodel.ub = fullmodel.ub';
    fullmodel.c = fullmodel.c';
    fullmodel.rxns = fullmodel.rxns';
    model.rev = model.rev';
    model.lb = model.lb';
    model.ub = model.ub';
    model.c = model.c';
    model.rxns = model.rxns';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Support function for expa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = clearrows(M)
    %Step 4: For all of the rows added to Tx in steps 2 and 3 check to make
    %sure that no row exists that is a non-negative combination of any other set
    %of rows in Tx. One method used is as follows:
    %Let A(i) equal the set of column indices, j, for which the elements of row
    %i equal zero. Then check to determine if there exists another row h for
    %which A(i) is a subset of A(h).
    %   A(i) is contained in A(h), i~=h
    %   where
    %   A(i) = {j:T(i,j)=0, 1<=j<=columns of T}
    %Thus if these equations hold true for any distinct os i and h, then row i
    %must be eliminated from T1
    %Here we will calculate the matrix test as follows:
    %In this matrix, if a column index is zero in row i and not zero in row h,
    %which would mean that A(i) is not contained in A(h), test(i,h) >= 1. So if
    %test(i,h) == 0 and i~=h, row i must be removed.
    %reverse test row. If test(i,h) == 1 remove row i
    test = double(M==0)*double(M'~=0);
    test = double(test==0) - eye(size(test,1));

    %Now we see if there are no rows with the exact same row indexes
    test = (test - tril(test)')>0;

    M = [M sum(test,2)];
    M = sortrows(M,size(M,2));
    x = length(find(sum(test,2)));
    M = M(1:size(M,1)-x,1:size(M,2)-1);
end

function model = deleteSrow(model,n)
    model.S(n,:) = [];
    if isfield(model,'mets')
        model.mets(n) = [];
    end
    if isfield(model,'metNames')
        model.metNames(n) = [];
    end
    if isfield(model,'metFormulas')
        model.metFormulas(n) = [];
    end
    if isfield(model,'metChEBIID')
        model.metChEBIID(n) = [];
    end
    if isfield(model,'metKEGGID')
        model.metKEGGID(n) = [];
    end
    if isfield(model,'metPubChemID')
        model.metPubChemID(n) = [];
    end
    if isfield(model,'metInChIString')
        model.metInChIString(n) = [];
    end
    if isfield(model,'metCharge')
        model.metCharge(n) = [];
    end
    if isfield(model,'b')
        model.b(n) = [];
    end
end
function model = deleteScol(model,n)
    model.S(:,n) = [];
    if isfield(model,'rxns')
        model.rxns(n) = [];
    end
    if isfield(model,'rev')
        model.rev(n) = [];
    end
    if isfield(model,'ub')
        model.ub(n) = [];
    end
    if isfield(model,'lb')
        model.lb(n) = [];
    end
    if isfield(model,'c')
        model.c(n) = [];
    end
    if isfield(model,'rules')
        model.rules(n) = [];
    end
    if isfield(model,'rxnGeneMat')
        model.rxnGeneMat(n,:) = [];
    end
    
    if isfield(model,'grRules')
        model.grRules(n) = [];
    end
    
    if isfield(model,'subSystems')
        model.subSystems(n) = [];
    end
    
    if isfield(model,'confidenceScores')
        model.confidenceScores(n) = [];
    end
    
    if isfield(model,'rxnReferences')
        model.rxnReferences(n) = [];
    end
    
    if isfield(model,'rxnECNumbers')
        model.rxnECNumbers(n) = [];
    end
    
    if isfield(model,'rxnNotes')
        model.rxnNotes(n) = [];
    end
    
    if isfield(model,'rxnNames')
        model.rxnNames(n) = [];
    end
end
function model = expa_elementary_deleteone(rx,model,original)
    %See if reaction is actually in the model
    if ~ismember(rx,model.rxns)
        fprintf('WARNING: Reaction not in the model\n')
        while ~strcmp(rx,'q') && ~ismember(rx,model.rxns)
            rx = input(['Input reaction again: (q to quit, p to print exchange' ...
                'reactions)\n'],'s');
            if strcmp(rx,'q')
                return
            elseif strcmp(rx,'p')
                x = find(findExcRxns(model));
                for i = 1:length(x)
                    fprintf([model.rxns{x(i)} '\t' num2str(model.lb(x(i))) '\t' ...
                       num2str(model.ub(x(i))) '\n'])
                end
            end
        end
    end
    %Create temporaray matrix without exchange reactions or metabolite in
    %question
    temp = full(model.S);
    y = find(strcmp(rx,model.rxns));
    x = find(temp(:,y));
    %find exchange reactions
    t = find(findExcRxns(model));
    %Delete from matrix exchange reactions and metabolite in question
    temp(:,t) = [];
    temp(x,:) = [];
    temp_rxns = model.rxns;
    temp_rxns(t) = [];
    %If by removeing the metabolite any other reaction becomes a sink or a
    %source (has only positive or negative stoichiometric coefficients), we
    %might want to delete those
    if ~isempty(find(sum(temp > 0)==0,1)) || ...
      ~isempty(find(sum(temp < 0)==0,1))
        fprintf(['WARNING: Reaction ' rx ' cannot simply be deleted\n'])
        fprintf('Reactions that would become sink or source are:\n')
        t = find(double(sum((temp>0))==0) + double(sum((temp<0))==0));
        printRxnFormula(model,temp_rxns(t));
        fprintf('Reactions were originally:\n')
        printRxnFormula(original,temp_rxns(t));

        %Get input on whether to delete them or not
        bol = askinput(['Would you like to delete these reactions?'...
            ' (y/n)'],{'y' 'n'});
        %If yes, delete them
        if bol == 'y'
            model = removeRxns(model,temp_rxns(t));
            %If no, return without deleting the exchange reaction
        else
            fprintf([rx ' was not deleted\n'])
            return
        end
    end
    %Delete exchange reaction and metabolite
    y = find(model.S(:,strcmp(rx,model.rxns)));
    model = removeMetabolites(model,model.mets{y});
    expa_elementary_analyze(model)
end
function model = expa_elementary_deleteone_met(met,model,original)
    %See if reaction is actually in the model
    if ~ismember(met,model.mets)
        fprintf('WARNING: Metabolite not in the model\n')
        while ~strcmp(met,'q') && ~ismember(met,model.mets)
            met = input('Input metabolite again: (q to quit)\n','s');
            if strcmp(met,'q')
                return
            end
        end
    end
    %Create temporaray matrix without metaboliet
    temp = full(model.S);
    y = find(strcmp(met,model.mets));
    %find exchange reactions
    t = find(findExcRxns(model));
    %Delete from matrix exchange reactions and metabolite in question
    temp(:,t) = [];
    temp(y,:) = [];
    temp_rxns = model.rxns;
    temp_rxns(t) = [];
    %If by removing the metabolite any other reaction becomes a sink or a
    %source (has only positive or negative stoichiometric coefficients), we
    %might want to delete those
    if ~isempty(find(sum(temp > 0)==0,1)) || ...
      ~isempty(find(sum(temp < 0)==0,1))
        fprintf(['WARNING: metabolite ' met ' cannot simply be deleted\n'])
        fprintf('Reactions that would become sink or source are:\n')
        t = find(double(sum((temp>0))==0) + double(sum((temp<0))==0));
        printRxnFormula(model,temp_rxns(t));
        fprintf('Reactions were originally:\n')
        printRxnFormula(original,temp_rxns(t));
        %Get input on whether to delete them or not
        bol = askinput(['Would you like to delete these reactions?'...
            ' (y/n)'],{'y' 'n'});
        %If yes, delete them
        if bol == 'y'
            model = removeRxns(model,temp_rxns(t));
            %If no, return without deleting the exchange reaction
        else
            fprintf([met ' was not deleted\n'])
            return
        end
    end
    %Delete exchange reaction and metabolite
    model = removeMetabolites(model,met);
    expa_elementary_analyze(model)
end





