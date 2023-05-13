Needs["xTools`xTension`", "xTension.wl"];
$DefInfoQ = False;

DefConstantSymbol[dimx];
DefManifold[MR, 1, {MR`r1, MR`r2, MR`r3, MR`r4}];
DefManifold[MX, dimx, IndexRangeNS[MX`a, MX`f]];
DefManifold[MF, {MR, MX}, IndexRangeNS[MF`A, MF`F]];
DefParameter[r];
DefConstantSymbol[L0];
DefConstantSymbol[Lt];
DefConstantSymbol[Lambda0];
DefMetric[1, metricMX[-MX`a, -MX`b], CDMX, OtherDependencies -> {r}, PrintAs -> "g", SymbolOfCovD -> {";", "\[Del]"}];
DefMetric[1, metricMF[-MF`A, -MF`B], CDMF, OtherDependencies -> {r}, PrintAs -> "g", SymbolOfCovD -> {";", "D"}];
DefTensor[V0[MF`A], MF];
DefTensor[T1[MX`a, -MX`b], {MX, r}];

DefCoordinateParameter[TangentMF -> TangentMR, r, eMR, edMR];

MUnit`BeginTestSection["FG expansion"];

If[!xTools`$SkipLongTest,
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
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["DefRiemannVarD"];

DefRiemannVarD[CDMF];
VerificationTest[
    VarD[RiemannCDMF[-MF`A, -MF`B, -MF`C, -MF`D]][Scalar[RicciCDMF[MF`A, MF`B] RicciCDMF[-MF`A, -MF`B]]]
    - 1/2 (
        metricMF[MF`B, MF`D] RicciCDMF[MF`A, MF`C]
        - metricMF[MF`B, MF`C] RicciCDMF[MF`A, MF`D]
        - metricMF[MF`A, MF`D] RicciCDMF[MF`B, MF`C]
        + metricMF[MF`A, MF`C] RicciCDMF[MF`B, MF`D]
    ) // ContractMetric // Simplification,
    0
];
VerificationTest[
    VarD[RiemannCDMF[-MF`A, -MF`B, -MF`C, -MF`D]][RicciScalarCDMF[]]
    - 1/2 (- metricMF[MF`A, MF`D] metricMF[MF`B, MF`C] + metricMF[MF`A, MF`C] metricMF[MF`B, MF`D]) // ContractMetric // Simplification,
    0
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["ETensor"];

VerificationTest[
    ETensor[metricMF[-MF`A, -MF`B]][MF`A, -MF`A]
,
    dimx + 1
];

VerificationTest[
    ETensor[metricMF[-MF`A, -MF`B]][MF`B, -MF`C] - delta[MF`B, -MF`C] // ToCanonical
,
    0
];

VerificationTest[
    ETensor[metricMF[-MF`A, -MF`B]] - ETensor[metricMF[-MF`C, -MF`D]]
,
    ETensor[0, {-MF`A, -MF`B}]
];

VerificationTest[
    ETensorProduct[ETensor[metricMF[-MF`A, -MF`B]], ETensor[metricMF[-MF`A, MF`B]]]
    - ETensor[metricMF[-MF`A, -MF`B] delta[-MF`C, MF`D], {-MF`A, -MF`B, -MF`C, MF`D}] // #[-MF`A, -MF`B, -MF`C, MF`D] &
,
    0
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["ETensor product with Scalar inside"];

VerificationTest[
    ETensorProduct[
        ETensor[Scalar[V0[MF`A] V0[-MF`A]] V0[MF`A], {MF`A}],
        ETensor[Scalar[V0[MF`A] V0[-MF`A]] V0[-MF`A], {-MF`A}]
    ] - ETensor[Scalar[V0[MF`A] V0[-MF`A]]^2 V0[MF`A] V0[-MF`B], {MF`A, -MF`B}] // ZeroETensorQ
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["ParamD and CovD"];

VerificationTest[
    SortCommParamDLeviCivitaCovD[ParamD[r]@CDMX[MX`a]@T1[MX`b, -MX`c]] - (
        ChristoffelToGradMetric[ChangeCovD[ParamD[r][metricMX[MX`a, MX`d] CDMX[-MX`d]@T1[MX`b, -MX`c]], CDMX, PD], metricMX] /. ParamD[p_]@PD[b_]@expr_ :> PD[b]@ParamD[p]@expr // ChangeCovD[#, PD, CDMX] &
    ) // ToCanonicalN // ContractMetric[#, AllowUpperDerivatives -> True] & // ToCanonicalN
,
    0
];

VerificationTest[
    ExpandParamDLeviCivitaChristoffel[ParamD[r]@ChristoffelCDMX[MX`a, -MX`b, -MX`c]] - (
        ChristoffelToGradMetric[ParamD[r]@ChristoffelCDMX[MX`a, -MX`b, -MX`c], metricMX] /. ParamD[p_]@PD[b_]@expr_ :> PD[b]@ParamD[p]@expr // ChangeCovD[#, PD, CDMX] &
    ) // ToCanonicalN
,
    0
];

VerificationTest[
    PdSymChristoffelToRiemann[ChangeCurvature[RiemannCDMF[-MF`A, -MF`B, -MF`C, MF`D], CDMF, PD]] - RiemannCDMF[-MF`A, -MF`B, -MF`C, MF`D] // ToCanonicalN
,
    0
];

MUnit`EndTestSection[];

UndefTensor /@ {eMR, edMR, V0, T1};
Undef /@ VisitorsOf@metricMX;
Undef /@ VisitorsOf@metricMF;
UndefMetric[metricMX];
UndefMetric[metricMF];
Undef /@ VisitorsOf@MX;
Undef /@ VisitorsOf@MR;
Undef /@ VisitorsOf@MF;
UndefManifold /@ {MF, MR, MX};
UndefConstantSymbol /@ {L0, Lt, Lambda0, dimx};
UndefParameter[r];