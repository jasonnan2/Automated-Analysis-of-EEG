function behModels = extractp(behModels,template,showall)
if nargin<3
    showall=1;
end
for i=1:height(behModels)
    model=behModels.model{i};
    if sum(strcmp(template,model.CoefficientNames)) & ~isempty(template)
        varName=template;
    else
        varName=append(template,behModels{i,1}{:});
    end
    neuralidx = find(strcmp(varName,model.CoefficientNames));
    pvals(i) = model.Coefficients.pValue(neuralidx);
end
behModels.(append(template,'neuralP'))=pvals';
if ~showall
    behModels(behModels.(append(template,'neuralP'))>0.05,:)=[];
end
end