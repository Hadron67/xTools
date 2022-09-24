<< xTension`;

MUnit`BeginTestSection["FG expansion"];

DefConstantSymbol[dimx];
DefManifold[MR, 1, {MR`r1, MR`r2, MR`r3, MR`r4}];
DefManifold[MX, dimx, {MX`a, MX`b, MX`c, MX`d, MX`e, MX`f}];
DefManifold[MF, {MR, MX}, {MF`A, MF`B, MF`C, MF`D, MF`E, MF`F}];
DefParameter[r];
DefConstantSymbol[L0];
DefConstantSymbol[Lt];
DefConstantSymbol[Lambda0];
DefMetric[1, metricMX[-MX`a, -MX`b], CDMX, OtherDependencies -> {r}, PrintAs -> "g", SymbolOfCovD -> {";", "\[Del]"}];
DefMetric[1, metricMF[-MF`A, -MF`B], CDMF, OtherDependencies -> {r}, PrintAs -> "g", SymbolOfCovD -> {";", "D"}];

DefCoordinateParameter[TangentMF -> TangentMR, r, eMR, edMR];
SetDecomposedRules[metricMF[-a_, -b_], DiagonalMatrix@{Lt^2/(4 r^2) edMR[-a] edMR[-b], 1/r metricMX[-a, -b]}];
SetDecomposedRules[metricMF[a_, b_], DiagonalMatrix@{(4 r^2)/Lt^2 eMR[a] eMR[b], r metricMX[a, b]}];

einsteinEom = (RicciCDMF[-MF`A, -MF`B] + (2 Lambda0)/(dimx - 1) metricMF[-MF`A, -MF`B]
    /. Lambda0 -> (dimx (dimx - 1))/(2 Lt^2)
    // ToDecomposed[CDMF]
    // SimplificationN
    // DropCoordinateBasis);

einsteinEomExpected = {{(metricMX[MX`b$16537$16669, MX`b$16538$16669]*(metricMX[MX`b$16539$16669, MX`b$16540$16669]*ParamD[r][metricMX[-MX`b$16537$16669, -MX`b$16539$16669]]*ParamD[r][metricMX[-MX`b$16538$16669, -MX`b$16540$16669]] -2*ParamD[r, r][metricMX[-MX`b$16537$16669, -MX`b$16538$16669]]))/4,(metricMX[MX`b$16548$16670, MX`b$16549$16670]*(-CDMX[-MX`b][ParamD[r][metricMX[-MX`b$16548$16670, -MX`b$16549$16670]]] +CDMX[-MX`b$16549$16670][ParamD[r][metricMX[-MX`b, -MX`b$16548$16670]]]))/2},{(metricMX[MX`b$16565$16671, MX`b$16566$16671]*(-CDMX[-MX`a][ParamD[r][metricMX[-MX`b$16565$16671, -MX`b$16566$16671]]] +CDMX[-MX`b$16566$16671][ParamD[r][metricMX[-MX`a, -MX`b$16565$16671]]]))/2,(Lt^2*RicciCDMX[-MX`a, -MX`b] + metricMX[MX`b$16582$16672, MX`b$16583$16672]*(2*r*ParamD[r][metricMX[-MX`a, -MX`b$16582$16672]]*ParamD[r][metricMX[-MX`b, -MX`b$16583$16672]] + metricMX[-MX`a, -MX`b]*ParamD[r][metricMX[-MX`b$16582$16672, -MX`b$16583$16672]]) + ParamD[r][metricMX[-MX`a, -MX`b]]*(-2 + dimx - r*metricMX[MX`b$16582$16672, MX`b$16583$16672]*ParamD[r][metricMX[-MX`b$16582$16672, -MX`b$16583$16672]]) -2*r*ParamD[r, r][metricMX[-MX`a, -MX`b]])/Lt^2}};

VerificationTest[einsteinEom - einsteinEomExpected // ContractMetric[#, metricMX] & // SimplificationN, {{0, 0}, {0, 0}}];

MUnit`EndTestSection[];

MUnit`BeginTestSection["ReplaceIndicesRules"];

DefManifold[M4, 4, {M4`a, M4`b, M4`c, M4`d, M4`e, M4`f}];
DefMetric[1, metricM4[-M4`a, -M4`b], CDM4, PrintAs -> "g"];

VerificationTest[ReplaceIndicesRules[RicciCDMF[-MF`A, -MF`B], TangentMF, TangentM4], {MF`A -> M4`a, MF`B -> M4`b}];

MUint`EndTestSection[];