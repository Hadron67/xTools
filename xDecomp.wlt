Needs["xTools`xDecomp`", "xDecomp.wl"];

$DefInfoQ = False;

DefConstantSymbol[dimx];
DefParameter /@ {r, t};
DefScalarFunction /@ {h, f};
DefManifold[Mf, dimx + 2, IndexRangeNS[Mf`A, Mf`H]];
DefManifold[Mx, dimx, IndexRangeNS[Mx`a, Mx`h]];
DefMetric[NoSignDet, metricMf[-Mf`A, -Mf`B], CDMf];
DefMetric[NoSignDet, metricMx[-Mx`a, -Mx`b], CDMx];

DefGBasis[decomp, {{et[Mf`A], edt[-Mf`A], t}, {er[Mf`A], edr[-Mf`A], r}}, {eX[Mx`a, -Mf`A]}];
DefTensor[n0[Mf`A], {Mf, r}];
DefTensor[v0[Mx`a], {Mx, r}];

DefConstantSymbol[L0];
DefParameter[r2];
DefManifold[Mr, dimx, IndexRangeNS[Mr`a, Mr`h]];
DefMetric[NoSignDet, metricMr[-Mr`a, -Mr`b], CDMr, OtherDependencies -> {r2}];
DefGBasis[fgc, {{er2[Mf`A], edr2[-Mf`A], r2}}, {eR[Mr`a, -Mf`A]}];

MUnit`BeginTestSection["Decomposing SSS metric ansatze"];

metricVal = CreateGCTensor[{
    {-t, -t} -> -h[r],
    {-r, -r} -> 1 / f[r],
    {-Mx`a, -Mx`b} -> r^2 metricMx[-Mx`a, -Mx`b]
}, -{decomp, decomp}];

metricInvVal = CreateGCTensor[{
    {t, t} -> -1 / h[r],
    {r, r} -> f[r],
    {Mx`a, Mx`b} -> 1 / r^2 metricMx[Mx`a, Mx`b]
}, {decomp, decomp}];
SetHeldMetric[covd1, metricMf, metricVal, metricInvVal];
AddCurvatureTensorsToHolder[covd1, decomp, ChristoffelCDMf];

VerificationTest[
    CachedGCTensor[covd1, RicciScalarCDMf, {}][],
    (1/(2 r^2 h[r]^2))(h[r] (2 h[
       r] (RicciScalarCDMx[] - dimx r f'[r]) -
     r^2 f'[r] h'[r]) +
  f[r] (-2 (-1 + dimx) dimx h[r]^2 + r^2 h'[r]^2 -
     2 r h[r] (dimx h'[r] + r h''[r])))
];

VerificationTest[
    (RiemannCDMf[-Mf`A, Mf`C, -Mf`B, Mf`D]
    * RiemannCDMf[-Mf`C, Mf`E, -Mf`D, Mf`F]
    * RiemannCDMf[-Mf`E, Mf`A, -Mf`F, Mf`B]
    /. GetAllHeldTensorRules[covd1]
    // ContractGCTensors[covd1]
    // NoScalar
    // ContractMetric
    // Simplification)
    - (1/(8 r^6 h[r]^3) (2 h[r]^3 (4 RiemannCDMx[-Mx`a, Mx`e, -Mx`c, Mx`f]RiemannCDMx[Mx`a, Mx`b, Mx`c, Mx`d]
       RiemannCDMx[-Mx`b, -Mx`e, -Mx`d, -Mx`f] +
      3 r^2 RicciScalarCDMx[] f'[r]^2) -
   2 (-1 + dimx) dimx f[r]^3 h[r] (4 (-2 + dimx) h[r]^2 + 3 r^2 h'[r]^2) - 
   3 f[r] (2 h[r]^3 (4  RicciCDMx[-Mx`a, -Mx`b]RicciCDMx[Mx`a, Mx`b] - 
         4  RiemannCDMx[-Mx`a, -Mx`c, -Mx`b, -Mx`d] RiemannCDMx[Mx`a, Mx`b, Mx`c, Mx`d] + (-1 + dimx) dimx r^2 f'[r]^2) + 
      dimx r^4 h[r] f'[r]^2 h'[r]^2) + 
   3 f[r]^2 (8 (-2 + dimx) h[r]^3 RicciScalarCDMx[] + 
      dimx r^4 f'[r] h'[r]^3 + 
      2 r^2 h[r] h'[
        r] (RicciScalarCDMx[] h'[r] - 
         dimx r^2 f'[r] h''[r])))
    ) // Simplification
,
    0
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["Outer product"];

covd1 /: HeldGCTensor[covd1, n0] = CreateGCTensor[{{r} -> Sqrt[f[r]]}, {decomp}];
VerificationTest[
    ETensor[(n0[Mf`A] n0[Mf`B] /. GetAllHeldTensorRules[covd1] // ContractGCTensors[covd1]) - CreateGCTensor[{{r, r} -> f[r]}, {decomp, decomp}][Mf`A, Mf`B], {Mf`A, Mf`B}] // ScreenDollarIndices
,
    CreateGCTensor[{}, {decomp, decomp}]
];

VerificationTest[
    n0[Mf`A] n0[Mf`B] n0[-Mf`A] n0[-Mf`B] /. GetAllHeldTensorRules[covd1] // ContractGCTensors[covd1]
,
    1
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["FG expansion"];

SetHeldMetric[fgcHolder, metricMf,
    CreateGCTensor[{
        {-r2, -r2} -> L0^2/(4 r2^2),
        {-Mr`a, -Mr`b} -> 1/r2 metricMr[-Mr`a, -Mr`b]
    }, -{fgc, fgc}],
    CreateGCTensor[{
        {r2, r2} -> (4 r2^2)/L0^2,
        {Mr`a, Mr`b} -> r2 metricMr[Mr`a, Mr`b]
    }, {fgc, fgc}]
];
AddCurvatureTensorsToHolder[fgcHolder, fgc, ChristoffelCDMf];
GCTensorHolderDAUseMetricVB[fgcHolder] ^= None;

fgEEom = ETensor[
    RicciCDMf[-Mf`A, -Mf`B] + dimx / L0^2 metricMf[-Mf`A, -Mf`B]
    /. GetAllHeldTensorRules[fgcHolder]
    // ContractGCTensors[fgcHolder], {-Mf`A, -Mf`B}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // Simplify;
fgEEomExpected = GCTensor[
    {{1/4  metricMr[Mr`a, Mr`b] (metricMr[Mr`c, Mr`d]  ParamD[r2][metricMr[-Mr`a, -Mr`c]] ParamD[r2][metricMr[-Mr`b, -Mr`d]] - 2 ParamD[r2, r2][metricMr[-Mr`a, -Mr`b]]), 
    ETensor[1/2  metricMr[Mr`b, Mr`c]  (-CDMr[-Mr`a][
    ParamD[r2][metricMr[-Mr`b, -Mr`c]]] + CDMr[-Mr`c][
    ParamD[r2][metricMr[-Mr`a, -Mr`b]]]), {-Mr`a}]}, 
    {ETensor[1/2 metricMr[Mr`b, Mr`c]  (-CDMr[-Mr`a][
    ParamD[r2][
    metricMr[-Mr`b, -Mr`c]]] + CDMr[-Mr`c][
    ParamD[r2][
    metricMr[-Mr`a, -Mr`b]]]), {-Mr`a}], 
    ETensor[(1/(
        L0^2))(L0^2  RicciCDMr[-Mr`a, -Mr`b] + 
        2 r2  metricMr[Mr`c, Mr`d]  ParamD[r2][
    metricMr[-Mr`a, -Mr`c]] ParamD[r2][
    metricMr[-Mr`b, -Mr`d]] + 
        metricMr[-Mr`a, -Mr`b]   metricMr[Mr`c, Mr`d] 
            ParamD[r2][
    metricMr[-Mr`c, -Mr`d]] + ParamD[r2][
    metricMr[-Mr`a, -Mr`b]] (-2 + 
            dimx - r2  metricMr[Mr`c, Mr`d] 
            ParamD[r2][
    metricMr[-Mr`c, -Mr`d]]) - 
        2 r2 ParamD[r2, r2][
    metricMr[-Mr`a, -Mr`b]]), {-Mr`a, -Mr`b}]}}
, {-fgc, -fgc}];

VerificationTest[
    Simplify@ToCanonical[fgEEom - fgEEomExpected, UseMetricOnVBundle -> None]
,
    CreateGCTensor[{}, -{fgc, fgc}]
]

MUnit`EndTestSection[];

MUnit`BeginTestSection["ZeroGCTensorQ"];

VerificationTest[ZeroGCTensorQ@CreateGCTensor[{}, -{fgc, fgc, fgc}]];

MUnit`EndTestSection[];

MUnit`BeginTestSection["PostETensorContract"];

VerificationTest[
    GCTensorContractTwo[CreateGCTensor[{
        {-Mx`a, -Mx`b} -> metricMx[-Mx`a, -Mx`b] + metricMx[-Mx`b, -Mx`a]
    }, -{decomp, decomp}], CreateGCTensor[{
        {Mx`a, Mx`b} -> metricMx[Mx`a, Mx`b] + metricMx[Mx`b, Mx`a]
    }, {decomp, decomp}], {1}, {1}, PostETensorContract -> ContractMetric]
    - CreateGCTensor[{{-Mx`a, Mx`b} -> 4 delta[-Mx`a, Mx`b]}, {-decomp, decomp}] // ZeroGCTensorQ
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["Covariant derivatives"];

v1 = CreateGCTensor[{
    {t} -> h[t, r],
    {r} -> f[t, r],
    {Mx`a} -> v0[Mx`a]
}, {decomp}];

t1 = CreateGCTensor[{
    {t, -t} -> h[t, r],
    {r, -r} -> f[t, r],
    {Mx`a, -Mx`b} -> v0[Mx`a] v0[-Mx`b]
}, {decomp, -decomp}];

chris = CachedGCTensor[covd1, ChristoffelCDMf, {1, -1, -1}];

VerificationTest[
    GCTensorPDDiv[v1, 1][] // NoScalar // ScreenDollarIndices
,
    D[f[t, r], r] + D[h[t, r], t] + PD[-Mx`a][v0[Mx`a]]
];

VerificationTest[
    0 == GCTensorCovDDiv[v1, 1, -decomp, chris][] - ContractGCTensors[
        PDGChart[-decomp, -Mf`A][v1[Mf`A]] + chris[Mf`B, -Mf`B, -Mf`A] v1[Mf`A]
    ,
        covd1
    ] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] &
];

VerificationTest[
    GCTensorCovDDiv[t1, 1, -decomp, chris] - ETensor[ContractGCTensors[
        PDGChart[-decomp, -Mf`B][t1[Mf`B, -Mf`A]]
        + chris[Mf`B, -Mf`B, -Mf`C] t1[Mf`C, -Mf`A]
        - chris[Mf`C, -Mf`B, -Mf`A] t1[Mf`B, -Mf`C]
    ,
        covd1
    ], {-Mf`A}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // ZeroGCTensorQ
];

VerificationTest[
    GCTensorCovDGrad[t1, -decomp, chris] - ETensor[ContractGCTensors[
        PDGChart[-decomp, -Mf`C][t1[Mf`A, -Mf`B]]
        + chris[Mf`A, -Mf`C, -Mf`D] t1[Mf`D, -Mf`B]
        - chris[Mf`D, -Mf`C, -Mf`B] t1[Mf`A, -Mf`D]
    ,
        covd1
    ], {Mf`A, -Mf`B, -Mf`C}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // ZeroGCTensorQ
];

MUnit`EndTestSection[];

UndefGBasis@decomp;
UndefGBasis@fgc;
UndefTensor /@ {n0, v0};
Undef /@ VisitorsOf@metricMf;
Undef /@ VisitorsOf@metricMx;
Undef /@ VisitorsOf@metricMr;
UndefMetric /@ {metricMx, metricMf, metricMr};
Undef /@ VisitorsOf@Mx;
Undef /@ VisitorsOf@Mf;
Undef /@ VisitorsOf@Mr;
UndefManifold /@ {Mx, Mf, Mr};
UndefConstantSymbol[dimx];
UndefConstantSymbol[L0];
UndefParameter /@ {r, t};