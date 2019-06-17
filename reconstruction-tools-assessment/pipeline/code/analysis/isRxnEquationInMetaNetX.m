function rxnWasFound = isRxnEquationInMetaNetX(equation,MNX)

genericEquation = makeGenericReaction(equation);
rxnWasFound = reactionFormulaInModel(MNX, genericEquation,1,1);

end